#include <unistd.h>

#include <cmath>
#include <cstdio>
#include <ctime>

#include "mirco_contactpredictors_kokkos.h"
#include "mirco_contactstatus_kokkos.h"
#include "mirco_evaluate_kokkos.h"
#include "mirco_matrixsetup_kokkos.h"
#include "mirco_nonlinearsolver_kokkos.h"
#include "mirco_warmstart_kokkos.h"

namespace MIRCO
{
  void Evaluate(double& pressure, const double Delta, const double LateralLength,
      const double GridSize, const double Tolerance, const int MaxIteration,
      const double CompositeYoungs, const bool WarmStartingFlag,
      const double ElasticComplianceCorrection, const ViewMatrix_d topology_d, const double zmax,
      const ViewVector_d meshgrid_d, const bool PressureGreenFunFlag)
  {
    // Initialise the area vector and force vector. Each element containing the
    // area and force calculated at every iteration.
    std::vector<double> totalForceVector;
    std::vector<double> contactAreaVector;
    double w_el = 0.0;

    // Initialise number of iterations
    int k = 0;

    // Points in contact in the previous iteration (only needed for warmstart)
    ViewVectorInt_d activeSetf_d;

    // Contact force at (xvf,yvf) predicted in the previous iteration.
    ViewVector_d pf_d;

    // Initialise the error in force
    double ErrorForce = std::numeric_limits<double>::max();

    while (ErrorForce > Tolerance && k < MaxIteration)
    {
      // Indices of the points predicted to be in contact
      ViewVectorInt_d activeSet0_d;
      // Coordinates of the points predicted to be in contact
      ViewVector_d xv0_d, yv0_d;
      // Indentation value of the half space at the predicted points of contact
      ViewVector_d b0_d;

      // First predictor for contact set
      ContactSetPredictor(
          activeSet0_d, xv0_d, yv0_d, b0_d, zmax, Delta, w_el, topology_d, meshgrid_d);

      // Initial number of predicted contact nodes.
      const int n0 = activeSet0_d.extent(0);

      ViewVector_d p0_d;
      // Second predictor for contact set
      // p_d --> contact forces at (xvf,yvf) predicted in the previous iteration but
      // are a part of currect predicted contact set. p_d is calculated in the
      // Warmstart function to be used in the NNLS to accelerate the simulation.
      if (WarmStartingFlag && k > 0)
      {
        p0_d = Warmstart(activeSet0_d, activeSetf_d, pf_d);
      }
      else
      {
        p0_d = ViewVector_d("p0_d", n0, 0.0);
      }

      auto H_d = MatrixGeneration::SetupMatrix(
          xv0_d, yv0_d, GridSize, CompositeYoungs, n0, PressureGreenFunFlag);

      // Defined as (u - u(bar)) in (Bemporad & Paggi, 2015)
      // Gap between the point on the topology and the half space
      // ViewVector_d w_d;

      // use Nonlinear solver --> Non-Negative Least Squares (NNLS) as in
      // (Bemporad & Paggi, 2015)
      nonlinearSolve(pf_d, activeSetf_d, p0_d, activeSet0_d, H_d, b0_d);

      // Compute total contact force and contact area
      double totalForce;
      double contactArea;
      ComputeContactForceAndArea(
          totalForce, contactArea, pf_d, GridSize, LateralLength, PressureGreenFunFlag);
      totalForceVector.push_back(totalForce);
      contactAreaVector.push_back(contactArea);

      // Elastic correction, used in the next iteration
      w_el = totalForce / ElasticComplianceCorrection;

      // Compute error due to nonlinear correction
      if (k > 0)
      {
        ErrorForce = abs(totalForceVector[k] - totalForceVector[k - 1]) / totalForceVector[k];
      }

      ++k;
    }

    if (ErrorForce > Tolerance)
      std::runtime_error("The solution did not converge in the maximum number of iterations");

    // Calculate the final force value at the end of the iteration.
    const double finalForce = totalForceVector[k - 1];

    // Mean pressure
    pressure = finalForce / (LateralLength * LateralLength);
  }

}  // namespace MIRCO
