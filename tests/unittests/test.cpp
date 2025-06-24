#include <gtest/gtest.h>
#include <stdlib.h>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialSymDenseMatrix.hpp>
#include <vector>

#include "../../src/mirco_filesystem_utils.h"
#include "../../src/mirco_linearsolver.h"
#include "../../src/mirco_nonlinearsolver.h"
#include "../../src/mirco_topology.h"
#include "../../src/mirco_warmstart.h"
#include "nonlinear_solver_test.h"

TEST(linearsolver, solves)
{
  int systemsize = 2;

  // Build the matrix
  Teuchos::SerialSymDenseMatrix<int, double> topology;
  topology.shape(systemsize);
  for (int i = 0; i < systemsize; i++)
  {
    topology(i, i) = 2;
    for (int j = 0; j < i; j++)
    {
      topology(i, j) = 1;
      topology(j, i) = 1;
    }
  }

  // Build the vectors
  Teuchos::SerialDenseVector<int, double> vector_x;
  Teuchos::SerialDenseVector<int, double> vector_b;

  // Bring right hand side into the correct form
  vector_b.size(systemsize);

  // Build right hand side
  for (int i = 0; i < systemsize; i++)
  {
    vector_b(i) = 1;
  }

  // Call linear solver
  vector_x = MIRCO::LinearSolver::Solve(topology, vector_b);

  EXPECT_NEAR(vector_x(0), 0.333333333333333, 1e-06);
  EXPECT_NEAR(vector_x(1), 0.333333333333333, 1e-06);
}

TEST_F(NonlinearSolverTest, primalvariable)
{
  y_ = MIRCO::NonLinearSolver::Solve(matrix_, b_vector_, x_vector_, w_);

  EXPECT_NEAR(y_(0), 163213.374921086, 1e-06);
  EXPECT_NEAR(y_(1), 43877.9231473546, 1e-06);
  EXPECT_NEAR(y_(2), 163702.923578063, 1e-06);
  EXPECT_NEAR(y_(3), 55159.5440853170, 1e-06);
  EXPECT_NEAR(y_(4), 10542.1713862417, 1e-06);
  EXPECT_NEAR(y_(5), 53809.0897795325, 1e-06);
  EXPECT_NEAR(y_(6), 148773.412150208, 1e-06);
  EXPECT_NEAR(y_(7), 83711.5732276221, 1e-06);
  EXPECT_NEAR(y_(8), 149262.960807186, 1e-06);
}

TEST_F(NonlinearSolverTest, dualvariable)
{
  MIRCO::NonLinearSolver::Solve(matrix_, b_vector_, x_vector_, w_);

  EXPECT_NEAR(w_(0, 0), 0, 1e-06);
  EXPECT_NEAR(w_(1, 0), 0, 1e-06);
  EXPECT_NEAR(w_(2, 0), 0, 1e-06);
  EXPECT_NEAR(w_(3, 0), 0, 1e-06);
  EXPECT_NEAR(w_(4, 0), 0, 1e-06);
  EXPECT_NEAR(w_(5, 0), 0, 1e-06);
  EXPECT_NEAR(w_(6, 0), 0, 1e-06);
  EXPECT_NEAR(w_(7, 0), 0, 1e-06);
  EXPECT_NEAR(w_(8, 0), 0, 1e-06);
}

TEST(FilesystemUtils, createrelativepath)
{
  std::string targetfilename = "input.dat";
  std::string sourcefilename = "../inputfiles/sourceinput.json";
  MIRCO::UTILS::ChangeRelativePath(targetfilename, sourcefilename);
  EXPECT_EQ(targetfilename, "../inputfiles/input.dat");
}

TEST(FilesystemUtils, keepabsolutpath)
{
  std::string targetfilename = "/root_dir/home/user/Input/input.dat";
  std::string sourcefilename = "../inputfiles/sourceinput.json";
  MIRCO::UTILS::ChangeRelativePath(targetfilename, sourcefilename);
  EXPECT_EQ(targetfilename, "/root_dir/home/user/Input/input.dat");
}

TEST(readtopology, RMG)
{
  int Resolution = 2;
  float HurstExponent = 0.1;
  bool RandomSeedFlag = false;
  int RandomGeneratorSeed = 95;
  Teuchos::SerialDenseMatrix<int, double> outsurf;
  int N = pow(2, Resolution);
  outsurf.shape(N + 1, N + 1);
  double InitialTopologyStdDeviation = 20.0;

  MIRCO::Rmg surface(
      Resolution, InitialTopologyStdDeviation, HurstExponent, RandomSeedFlag, RandomGeneratorSeed);
  surface.GetSurface(outsurf);

  EXPECT_NEAR(outsurf(0, 0), 23.5435469989256, 1e-06);
  EXPECT_NEAR(outsurf(0, 1), 30.2624522170979, 1e-06);
  EXPECT_NEAR(outsurf(0, 2), 69.5813622417479, 1e-06);
  EXPECT_NEAR(outsurf(0, 3), 43.5026425381265, 1e-06);
  EXPECT_NEAR(outsurf(0, 4), 23.5435469989256, 1e-06);
  EXPECT_NEAR(outsurf(1, 0), 68.8507553267314, 1e-06);
  EXPECT_NEAR(outsurf(1, 1), 73.8350740079714, 1e-06);
  EXPECT_NEAR(outsurf(1, 2), 77.9927972851754, 1e-06);
  EXPECT_NEAR(outsurf(1, 3), 35.2927793006724, 1e-06);
  EXPECT_NEAR(outsurf(1, 4), 22.6620325442329, 1e-06);
  EXPECT_NEAR(outsurf(2, 0), 39.1583562054882, 1e-06);
  EXPECT_NEAR(outsurf(2, 1), 19.2247183888878, 1e-06);
  EXPECT_NEAR(outsurf(2, 2), 79.1711886771701, 1e-06);
  EXPECT_NEAR(outsurf(2, 3), 5.66729306836534, 1e-06);
  EXPECT_NEAR(outsurf(2, 4), 41.3691438722521, 1e-06);
  EXPECT_NEAR(outsurf(3, 0), 59.1811726494348, 1e-06);
  EXPECT_NEAR(outsurf(3, 1), 21.2400598989696, 1e-06);
  EXPECT_NEAR(outsurf(3, 2), 54.6656122080671, 1e-06);
  EXPECT_NEAR(outsurf(3, 3), 28.0246974768169, 1e-06);
  EXPECT_NEAR(outsurf(3, 4), 6.72730409669533, 1e-06);
  EXPECT_NEAR(outsurf(4, 0), 23.5435469989256, 1e-06);
  EXPECT_NEAR(outsurf(4, 1), 0, 1e-03);
  EXPECT_NEAR(outsurf(4, 2), 30.6777944575233, 1e-06);
  EXPECT_NEAR(outsurf(4, 3), 35.2191824993355, 1e-06);
  EXPECT_NEAR(outsurf(4, 4), 23.5435469989256, 1e-06);
}

TEST(warmstarting, warmstart)
{
  Teuchos::SerialDenseMatrix<int, double> x0;
  std::vector<double> xv0, yv0, xvf, yvf, pf;

  xv0.resize(3);
  yv0.resize(3);
  x0.shape(3, 1);
  xvf.resize(2);
  yvf.resize(2);
  pf.resize(2);

  xv0[0] = 1;
  xv0[1] = 3;
  xv0[2] = 5;

  yv0[0] = 2;
  yv0[1] = 4;
  yv0[2] = 6;

  xvf[0] = 1;
  xvf[1] = 5;

  yvf[0] = 2;
  yvf[1] = 6;

  pf[0] = 10;
  pf[1] = 30;

  MIRCO::Warmstart(x0, xv0, yv0, xvf, yvf, pf);

  EXPECT_EQ(x0(0, 0), 10);
  EXPECT_EQ(x0(1, 0), 0);
  EXPECT_EQ(x0(2, 0), 30);
}

TEST(warmstarting, warmstart2)
{
  Teuchos::SerialDenseMatrix<int, double> x0;
  std::vector<double> xv0, yv0, xvf, yvf, pf;

  xv0.resize(3);
  yv0.resize(3);
  x0.shape(3, 1);
  xvf.resize(4);
  yvf.resize(4);
  pf.resize(4);

  xv0[0] = 1;
  xv0[1] = 3;
  xv0[2] = 5;

  yv0[0] = 2;
  yv0[1] = 4;
  yv0[2] = 6;

  xvf[0] = 1;
  xvf[1] = 7;
  xvf[2] = 9;
  xvf[3] = 5;


  yvf[0] = 2;
  yvf[1] = 8;
  yvf[2] = 10;
  yvf[3] = 6;

  pf[0] = 10;
  pf[1] = 50;
  pf[2] = 70;
  pf[3] = 30;

  MIRCO::Warmstart(x0, xv0, yv0, xvf, yvf, pf);

  EXPECT_EQ(x0(0, 0), 10);
  EXPECT_EQ(x0(1, 0), 0);
  EXPECT_EQ(x0(2, 0), 30);
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
