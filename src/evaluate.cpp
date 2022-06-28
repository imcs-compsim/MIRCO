#include <Epetra_SerialSpdDenseSolver.h>
#include <Epetra_SerialSymDenseMatrix.h>
#include <omp.h>
#include <unistd.h>
#include <Teuchos_Assert.hpp>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
using namespace std;
#include "contactpredictors.h"
#include "contactstatus.h"
#include "evaluate.h"
#include "matrixsetup.h"
#include "nonlinearsolver.h"


void MIRCO::Evaluate(double &force, double Delta, double lato, double delta, double errf,
    double to1, int max_iter, double E, bool flagwarm, double k_el,
    Epetra_SerialDenseMatrix topology, double zmax, std::vector<double> x)
{
  omp_set_num_threads(6);  // 6 seems to be optimal

  auto start = std::chrono::high_resolution_clock::now();

  time_t now = time(0);
  tm *ltm = localtime(&now);
  std::cout << "Time is: " << ltm->tm_hour << ":";
  std::cout << 1 + ltm->tm_min << endl;

  vector<double> area0;
  vector<double> force0;
  double w_el = 0, area = 0;
  int k = 0, n0 = 0;
  std::vector<double> xv0, yv0, b0, xvf, yvf, pf;
  Epetra_SerialDenseMatrix x0;
  int nf = 0;
  Epetra_SerialDenseMatrix A;
  Epetra_SerialDenseMatrix y;

  while (errf > to1 && k < max_iter)
  {
    // First predictor for contact set
    // All points, for which gap is bigger than the displacement of the rigid
    // indenter, cannot be in contact and thus are not checked in nonlinear solve
    // @{
    MIRCO::ContactSetPredictor(n0, xv0, yv0, b0, zmax, Delta, w_el, x, topology);

    A.Shape(xv0.size(), xv0.size());

    // Construction of the Matrix H = A
    MIRCO::MatrixGeneration matrix1;
    matrix1.SetUpMatrix(A, xv0, yv0, delta, E, n0);

    // Second predictor for contact set
    // @{
    MIRCO::InitialGuessPredictor(flagwarm, k, n0, nf, xv0, yv0, pf, x0, b0, xvf, yvf);
    // }

    // {
    Epetra_SerialDenseMatrix b0new;
    b0new.Shape(b0.size(), 1);
    // } Parallel region makes around this makes program slower

#pragma omp parallel for schedule(static, 16)  // Always same workload -> Static!
    for (long unsigned int i = 0; i < b0.size(); i++)
    {
      b0new(i, 0) = b0[i];
    }

    Epetra_SerialDenseMatrix w;

    MIRCO::NonLinearSolver solution2;
    solution2.NonlinearSolve(A, b0new, x0, w, y);  // y -> sol, w -> wsol; x0 -> y0

    // Compute number of contact node
    // @{
    MIRCO::ComputeContactNodes(xvf, yvf, pf, nf, y, xv0, yv0);
    // }

    // Compute contact force and contact area
    // @{
    MIRCO::ComputeContactForceAndArea(force0, area0, w_el, nf, pf, k, delta, lato, k_el);
    // }

    // Compute error due to nonlinear correction
    // @{
    if (k > 0)
    {
      errf = abs(force0[k] - force0[k - 1]) / force0[k];

      // errw(k) = abs((w_el0(k+1)-w_el0(k))/w_el0(k+1));
      // It appears that this is only a debugging variable without any uses,
      // therefore im not gonna implement this here.
    }
    k += 1;
    // }
  }

  TEUCHOS_TEST_FOR_EXCEPTION(errf > to1, std::out_of_range,
      "The solution did not converge in the maximum number of iternations defined");
  // @{

  force = force0[k - 1];
  area = area0[k - 1];

  // Mean pressure
  double sigmaz = force / pow(lato, 2);
  // Pressure unit per depth
  double pressz = sigmaz;
  cout << "k= " << k << " nf= " << nf << endl;
  cout << "Force= " << force << endl;
  cout << "area= " << area << endl;
  cout << "Mean pressure is:" + std::to_string(sigmaz) +
              " ; pressure unit per depth is:" + std::to_string(pressz) + " . \n";

  auto finish = std::chrono::high_resolution_clock::now();
  double elapsedTime2 = std::chrono::duration_cast<std::chrono::seconds>(finish - start).count();
  std::cout << "Elapsed time is: " + to_string(elapsedTime2) + "s." << endl;
}
