#include "mirco_evaluate.h"
#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <omp.h>
#include <unistd.h>
#include <Teuchos_Assert.hpp>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "mirco_contactpredictors.h"
#include "mirco_contactstatus.h"
#include "mirco_matrixsetup.h"
#include "mirco_nonlinearsolver.h"


void MIRCO::Evaluate(double& pressure, double Delta, double LateralLength, double GridSize,
    double Tolerance, int MaxIteration, double CompositeYoungs, bool WarmStartingFlag,
    double ElasticComplianceCorrection, Teuchos::SerialDenseMatrix<int,double>& topology, double zmax,
    std::vector<double>& meshgrid)
{
  omp_set_num_threads(6);  // 6 seems to be optimal

  // Initialise the area vector and force vector. Each element containing the
  // area and force calculated at every iteration.
  std::vector<double> area0;
  std::vector<double> force0;
  double w_el = 0.0;

  // Initialise number of iteration, k, and initial number of predicted contact
  // nodes, n0.
  int k = 0, n0 = 0;

  // Coordinates of the points predicted to be in contact.
  std::vector<double> xv0, yv0;
  // Coordinates of the points in contact in the previous iteration.
  std::vector<double> xvf, yvf;
  // Indentation value of the half space at the predicted points of contact.
  std::vector<double> b0;
  // Contact force at (xvf,yvf) predicted in the previous iteration.
  std::vector<double> pf;

  // x0 --> contact forces at (xvf,yvf) predicted in the previous iteration but
  // are a part of currect predicted contact set. x0 is calculated in the
  // Warmstart function to be used in the NNLS to accelerate the simulation.
  Teuchos::SerialDenseMatrix<int,double> x0;

  // The number of nodes in contact in the previous iteration.
  int nf = 0;

  // The influence coefficient matrix (Discrete version of Green Function)
  Teuchos::SerialDenseMatrix<int,double> A;
  // Solution containing force
  Teuchos::SerialDenseMatrix<int,double> y;

  // Initialise the error in force
  double ErrorForce = std::numeric_limits<double>::max();
  while (ErrorForce > Tolerance && k < MaxIteration)
  {
    // First predictor for contact set
    // @{
    MIRCO::ContactSetPredictor(n0, xv0, yv0, b0, zmax, Delta, w_el, meshgrid, topology);

    A.shape(xv0.size(), xv0.size());

    // Construction of the Matrix A
    MIRCO::MatrixGeneration matrix1;
    matrix1.SetUpMatrix(A, xv0, yv0, GridSize, CompositeYoungs, n0);

    // Second predictor for contact set
    // @{
    MIRCO::InitialGuessPredictor(WarmStartingFlag, k, n0, nf, xv0, yv0, pf, x0, b0, xvf, yvf);
    // }

    // {
    Teuchos::SerialDenseMatrix<int,double> b0new;
    b0new.shape(b0.size(), 1);
    // } Parallel region makes around this makes program slower
#pragma omp parallel for schedule(static, 16)  // Always same workload -> Static!
    for (long unsigned int i = 0; i < b0.size(); i++)
    {
      b0new(i, 0) = b0[i];
    }

    // Defined as (u - u(bar)) in (Bemporad & Paggi, 2015)
    // Gap between the point on the topology and the half space
    Teuchos::SerialDenseMatrix<int,double> w;

    // use Nonlinear solver --> Non-Negative Least Squares (NNLS) as in
    // (Bemporad & Paggi, 2015)
    MIRCO::NonLinearSolver solution2;
    solution2.NonlinearSolve(A, b0new, x0, w, y);

    // Compute number of contact node
    // @{
    MIRCO::ComputeContactNodes(xvf, yvf, pf, nf, y, xv0, yv0);
    // }

    // Compute contact force and contact area
    // @{
    MIRCO::ComputeContactForceAndArea(
        force0, area0, w_el, nf, pf, k, GridSize, LateralLength, ElasticComplianceCorrection);
    // }

    // Compute error due to nonlinear correction
    // @{
    if (k > 0)
    {
      ErrorForce = abs(force0[k] - force0[k - 1]) / force0[k];
    }
    k += 1;
    // }
  }

  TEUCHOS_TEST_FOR_EXCEPTION(ErrorForce > Tolerance, std::out_of_range,
      "The solution did not converge in the maximum number of iternations defined");
  // @{

  // Calculate the final force and area value at the end of the iteration.
  const double force = force0[k - 1];
  const double area = area0[k - 1];

  // Mean pressure
  double sigmaz = force / pow(LateralLength, 2);
  pressure = sigmaz;
}
