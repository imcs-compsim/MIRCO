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

#include "mirco_contactpredictors.h"
#include "mirco_contactstatus.h"
#include "mirco_matrixsetup_kokkos.h"
#include "mirco_nonlinearsolver.h"

void MIRCO::Evaluate(double& pressure, const double Delta, const double LateralLength,
    const double GridSize, const double Tolerance, const int MaxIteration,
    const double CompositeYoungs, const bool WarmStartingFlag,
    const double ElasticComplianceCorrection,
    const Teuchos::SerialDenseMatrix<int, double>& topology, const double zmax,
    const std::vector<double>& meshgrid, const bool PressureGreenFunFlag)
{
  // Initialise the area vector and force vector. Each element containing the
  // area and force calculated at every iteration.
  std::vector<double> area0;
  std::vector<double> force0;
  double w_el = 0.0;

  // Initialise number of iteration, k, and initial number of predicted contact
  // nodes, n0.
  int k = 0;  ///, n0 = 0;//#local to while scope

  // Coordinates of the points predicted to be in contact.
  /// std::vector<double> xv0, yv0;//# this should not be declared here as it is local to the while
  /// loop
  // Coordinates of the points in contact in the previous iteration.
  std::vector<double> xvf, yvf;
  // Indentation value of the half space at the predicted points of contact.
  std::vector<double> b0;
  // Contact force at (xvf,yvf) predicted in the previous iteration.
  std::vector<double> pf;

  // The number of nodes in contact in the previous iteration.
  int nf = 0;

  // The influence coefficient matrix (Discrete version of Green Function)
  ///  ViewMatrix_h H;//==A
  // Solution containing force
  ViewVector_h p_star;  //==y; //# <-- means previously called

  // Initialise the error in force
  double ErrorForce = std::numeric_limits<double>::max();
  while (ErrorForce > Tolerance && k < MaxIteration)
  {
    // Initial number of predicted contact nodes.
    int n0;
    // Coordinates of the points predicted to be in contact.
    std::vector<double> xv0, yv0;  // # change to ViewVector_d/h xv0, yv0;
    // // First predictor for contact set
    // # xv0 and xy0 and b0 are cleared, so are only out and not in; n0 is only out
    MIRCO::ContactSetPredictor(n0, xv0, yv0, b0, zmax, Delta, w_el, meshgrid, topology);

    // Second predictor for contact set
    // @{
    // x0 --> contact forces at (xvf,yvf) predicted in the previous iteration but
    // are a part of currect predicted contact set. x0 is calculated in the
    // Warmstart function to be used in the NNLS to accelerate the simulation.
    Teuchos::SerialDenseMatrix<int, double>
        x0;  // #ViewVector_h x0; //# this can be scoped in while, as it always gets reshaped in
             // InitialGuessPredictor() //# in fact we should just make InitialGuessPredictor()
             // return x0, as it is the only output

    MIRCO::InitialGuessPredictor(WarmStartingFlag, k, n0, xv0, yv0, pf, x0, b0, xvf, yvf);
    // }

    // Construction of the Matrix A
    /*if(!HOSTONLY) { //# or something like this?
      ViewMatrix_d H_d;
      Kokkos::deep_copy(H_d, H_h);
      H = H_d;
    }*/
    auto xv0_d = toKokkos(xv0);
    auto yv0_d = toKokkos(yv0);
    // ViewMatrix_d H_d("H", xv0.size(), xv0.size());
    auto H_d = MIRCO::MatrixGeneration::SetupMatrix(
        xv0_d, yv0_d, GridSize, CompositeYoungs, n0, PressureGreenFunFlag);

    /// ViewMatrix_h H_h = Kokkos::create_mirror_view_and_copy(MemorySpace_Host_t(), H_d); //# in
    /// future maybe

    Teuchos::SerialDenseMatrix A = kokkosMatrixToTeuchos(H_d);



    // {
    // Defined as (u - u(bar)) in (Bemporad & Paggi, 2015)
    // Gap between the point on the topology and the half space
    Teuchos::SerialDenseMatrix<int, double> w;

    // use Nonlinear solver --> Non-Negative Least Squares (NNLS) as in
    // (Bemporad & Paggi, 2015)
    auto y = MIRCO::NonLinearSolver::Solve(A, b0, x0, w);
    // #Q Why is x0 (or y0...) a matrix with 1 column and not just a vector. there is no reason i
    // think

    // Compute number of contact node
    MIRCO::ComputeContactNodes(xvf, yvf, pf, nf, y, xv0, yv0);

    // Compute contact force and contact area
    MIRCO::ComputeContactForceAndArea(force0, area0, w_el, nf, pf, k, GridSize, LateralLength,
        ElasticComplianceCorrection, PressureGreenFunFlag);

    // Compute error due to nonlinear correction
    if (k > 0)
    {
      ErrorForce = abs(force0[k] - force0[k - 1]) / force0[k];
    }
    k += 1;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(ErrorForce > Tolerance, std::out_of_range,
      "The solution did not converge in the maximum number of iternations defined");

  // Calculate the final force value at the end of the iteration.
  const double force = force0[k - 1];

  // Mean pressure
  double sigmaz = force / (LateralLength * LateralLength);
  pressure = sigmaz;
}
