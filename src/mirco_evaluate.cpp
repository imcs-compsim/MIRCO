#include "mirco_evaluate.h"

#include <unistd.h>

#include <cmath>
#include <cstdio>
#include <ctime>

#include "mirco_contactpredictors.h"
#include "mirco_contactstatus.h"
#include "mirco_matrixsetup.h"
#include "mirco_nonlinearsolver.h"
#include "mirco_warmstart.h"

namespace MIRCO
{
  void Evaluate(double& pressure, double& effectiveContactAreaFraction, const double Delta,
      const double LateralLength, const double GridSize, const double Tolerance,
      const int MaxIteration, const double CompositeYoungs, const bool WarmStartingFlag,
      const double ElasticComplianceCorrection, const ViewMatrix_d topology, const double zmax,
      const ViewVector_d meshgrid, const bool PressureGreenFunFlag)
  {
    // Initialise the area vector and force vector. Each element contains the
    // area and force calculated at every iteration.
    std::vector<double> totalForceVector;
    std::vector<double> contactAreaVector;
    double w_el = 0.0;

    // Initialise number of iterations
    int k = 0;

    // Points in contact in the previous iteration (only needed for warmstart)
    ViewVectorInt_d activeSetf;

    // Contact force at (xvf,yvf) predicted in the previous iteration
    ViewVector_d pf;

    // Difference in total force between current and previous iteration; used as a convergence
    // criterion
    double deltaTotalForce = std::numeric_limits<double>::max();

    while (deltaTotalForce > Tolerance && k < MaxIteration)
    {
      // Indices of the points predicted to be in contact
      ViewVectorInt_d activeSet0;
      // Coordinates of the points predicted to be in contact
      ViewVector_d xv0, yv0;
      // Indentation value of the half space at the predicted points of contact
      ViewVector_d b0;

      // First predictor for contact set
      ContactSetPredictor(activeSet0, xv0, yv0, b0, zmax, Delta, w_el, topology, meshgrid);

      // Initial number of predicted contact nodes.
      const int n0 = activeSet0.extent(0);

      ViewVector_d p0;
      if (WarmStartingFlag && k > 0)
      {
        // Warmstart
        // p0 --> contact forces at (xvf,yvf) predicted in the previous iteration but
        // are a part of currect predicted contact set. p0 is calculated in the
        // Warmstart function to be used in the NNLS to accelerate the simulation.
        p0 = Warmstart(activeSet0, activeSetf, pf);
      }
      else
      {
        p0 = ViewVector_d("p0_d", n0, 0.0);
      }

      auto H = MatrixGeneration::SetupMatrix(
          xv0, yv0, GridSize, CompositeYoungs, n0, PressureGreenFunFlag);

      // Defined as (u - u(bar)) in (Bemporad & Paggi, 2015)
      // Gap between the point on the topology and the half space
      // ViewVector_d w;

      // use Nonlinear solver --> Non-Negative Least Squares (NNLS) as in
      // (Bemporad & Paggi, 2015)
      nonlinearSolve(pf, activeSetf, p0, activeSet0, H, b0);

      // Compute total contact force and contact area
      double totalForce;
      double contactArea;
      ComputeContactForceAndArea(
          totalForce, contactArea, pf, GridSize, LateralLength, PressureGreenFunFlag);
      totalForceVector.push_back(totalForce);
      contactAreaVector.push_back(contactArea);

      // Elastic correction, used in the next iteration
      w_el = totalForce / ElasticComplianceCorrection;

      // Compute error due to nonlinear correction
      if (k > 0)
      {
        deltaTotalForce = abs(totalForceVector[k] - totalForceVector[k - 1]) / totalForceVector[k];
      }

      ++k;
    }

    if (deltaTotalForce > Tolerance)
      std::runtime_error("The solver did not converge in the maximum number of iterations.");

    // Calculate the final force value at the end of the iteration
    const double finalForce = totalForceVector.back();

    // Mean pressure
    pressure = finalForce / (LateralLength * LateralLength);

    // Effective contact area in converged state
    effectiveContactAreaFraction = contactAreaVector.back() / (LateralLength * LateralLength);
  }

}  // namespace MIRCO
