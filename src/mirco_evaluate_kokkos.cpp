#include "mirco_evaluate_kokkos.h"

#include <omp.h>
#include <unistd.h>

#include <Teuchos_Assert.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseSolver.hpp>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <vector>

#include "mirco_contactpredictors_kokkos.h"
#include "mirco_contactstatus_kokkos.h"
#include "mirco_matrixsetup_kokkos.h"
#include "mirco_nonlinearsolver_kokkos.h"
#include "mirco_warmstart_kokkos.h"

double MIRCO::Evaluate(const double Delta, const double LateralLength,
    const double GridSize, const double Tolerance, const int MaxIteration,
    const double CompositeYoungs, const bool WarmStartingFlag,
    const double ElasticComplianceCorrection,
    const Teuchos::SerialDenseMatrix<int, double>& topology, const double zmax,
    const std::vector<double>& meshgrid, const bool PressureGreenFunFlag)
{
  // Initialise the area vector and force vector. Each element containing the
  // area and force calculated at every iteration.
  std::vector<double> contactAreaVector;
  std::vector<double> totalForceVector;
  double w_el = 0.0;

  std::cout << "---------------------------Evaluate Kokkos\n";

  // Initialise number of iteration, k, and initial number of predicted contact
  // nodes, n0.
  int k = 0;  ///, n0 = 0;//#local to while scope

  // Coordinates of the points predicted to be in contact.
  /// std::vector<double> xv0, yv0;//# this should not be declared here as it is local to the while
  /// loop
  // Coordinates of the points in contact in the previous iteration.
  ViewVector_d xvf_d, yvf_d;
  // Contact force at (xvf,yvf) predicted in the previous iteration.
  ViewVector_d pf_d;

  // The number of nodes in contact in the previous iteration.
  int nf = 0;

  // The influence coefficient matrix (Discrete version of Green Function)
  ///  ViewMatrix_h H;//==A
  // Solution containing force
  /// ViewVector_h p_star;  //==y; //# <-- means previously called

  ViewVector_d meshgrid_d = toKokkos(meshgrid);
  ViewMatrix_d topology_d = toKokkos(topology);

  // Initialise the error in force
  double ErrorForce = std::numeric_limits<double>::max();  // # can just do tolerance + 1.0 so you
                                                           // dont need to include numeric_limits
  while (ErrorForce > Tolerance && k < MaxIteration)
  {
    // Initial number of predicted contact nodes.
    int n0;
    // Coordinates of the points predicted to be in contact.
    ViewVector_d xv0_d, yv0_d;

    // Indentation value of the half space at the predicted points of contact.
    // b0(i) = Delta + w_el - (zmax - topology(ri, ci));
    //                        ^\xi_{max}  ^ \xi             in Bemporad 2015
    // w_el = force0[k] / ElasticComplianceCorrection;
    // elastic_compliance_correction_ = LateralLength * composite_youngs_ / shape_factor_; where
    // LateralLength is the side length of the whole big (FEM) element (projected onto boundary),
    // which contains a bunch of smaller BEM elements (see Mayr 2019)
    // cf: \overbar{u} would be Delta - (zmax - topology(ri, ci))         without the w_el
    ViewVector_d b0_d;  // # change to ViewVector_d/h xv0, yv0;

    // // First predictor for contact set
    // # xv0 and xy0 and b0 are cleared, so are only out and not in; n0 is only out
    MIRCO::ContactSetPredictor(n0, xv0_d, yv0_d, b0_d, zmax, Delta, w_el, meshgrid_d, topology_d);

    // Second predictor for contact set
    // @{
    // x0 --> contact forces at (xvf,yvf) predicted in the previous iteration but
    // are a part of currect predicted contact set. x0 is calculated in the
    // Warmstart function to be used in the NNLS to accelerate the simulation.
    Teuchos::SerialDenseMatrix<int, double>
        x0;  // #ViewVector_h x0; //# this can be scoped in while, as it always gets reshaped in
             // InitialGuessPredictor() //# in fact we should just make InitialGuessPredictor()
             // return x0, as it is the only output
             
             
    ViewVector_d p0p_d;
    if (WarmStartingFlag && k > 0)
    {
      ViewVector_h p0_h = MIRCO::Warmstart(Kokkos::create_mirror_view_and_copy(MemorySpace_Host_t(), xv0_d),
        Kokkos::create_mirror_view_and_copy(MemorySpace_Host_t(), yv0_d),
        Kokkos::create_mirror_view_and_copy(MemorySpace_Host_t(), xvf_d),
        Kokkos::create_mirror_view_and_copy(MemorySpace_Host_t(), yvf_d),
        Kokkos::create_mirror_view_and_copy(MemorySpace_Host_t(), pf_d));
        
      p0p_d = Kokkos::create_mirror_view_and_copy(MemorySpace_ofDefaultExec_t(), p0_h);
    }
    else
    {
      if (b0_d.extent(0) > 0)
      {
        p0p_d = ViewVector_d("InitialGuessPredictor", n0, 0.0);
      }
    }

        
    // }
    auto H_d = MIRCO::MatrixGeneration::SetupMatrix(
        xv0_d, yv0_d, GridSize, CompositeYoungs, n0, PressureGreenFunFlag);


        //## we have at this point used the indexing of xv0 and such. index 0 of a quantity like p0_d is wherever that xv0 active set starts

    // {
    // Defined as (u - u(bar)) in (Bemporad & Paggi, 2015)
    // Gap between the point on the topology and the half space
    ///ViewVector_d w_d;
    
    int activeSetSize;
    // use Nonlinear solver --> Non-Negative Least Squares (NNLS) as in
    // (Bemporad & Paggi, 2015)
    //## here, we make an active set that holds in no particular order the references to the indices. but the thing is that p_d inside of this function, as the non-compacted vector, still holds
    MIRCO::NonLinearSolver::solve(H_d, b0_d, p0p_d, activeSetSize);///w_d);
    // #Q Why is x0 (or y0...) a matrix with 1 column and not just a vector. there is no reason i
    // think

    // Compute number of contact node
    //# TODO: this function, I think, we can actually just integrate into nonlinear solve and it will all be more efficient because we already have a compact form basically in nonlinear solve
    MIRCO::ComputeContactNodes(xvf_d, yvf_d, pf_d, activeSetSize, p0p_d, xv0_d, yv0_d);

    // Compute total contact force and contact area
    double totalForce;
    double contactArea;
    MIRCO::ComputeContactForceAndArea(totalForce, contactArea, pf_d, GridSize, LateralLength, PressureGreenFunFlag);
    
    totalForceVector.push_back(totalForce);
    contactAreaVector.push_back(contactArea);
    
    // Elastic correction, used in the next iteration
    w_el = totalForce / ElasticComplianceCorrection;

    // Compute error due to nonlinear correction
    if (k > 0)
    {
      ErrorForce = abs(totalForceVector[k] - totalForceVector[k - 1]) / totalForceVector[k];
    }
    k += 1;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(ErrorForce > Tolerance, std::out_of_range,
      "The solution did not converge in the maximum number of iternations defined");

  // Calculate the final force value at the end of the iteration.
  const double finalForce = totalForceVector[k - 1];

  // Mean pressure
  return finalForce / (LateralLength * LateralLength);
}
