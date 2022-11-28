#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialSymDenseMatrix.h>
#include <gtest/gtest.h>
#include <stdlib.h>
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
  Epetra_SerialSymDenseMatrix topology;
  topology.Shape(systemsize);
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
  Epetra_SerialDenseMatrix vector_x;
  Epetra_SerialDenseMatrix vector_b;

  // Bring matrices in correct form
  vector_x.Shape(systemsize, 1);
  vector_b.Shape(systemsize, 1);

  // Build right hand side
  for (int i = 0; i < systemsize; i++)
  {
    vector_b(i, 0) = 1;
  }

  // Call linear solver
  MIRCO::LinearSolver linearsolver;
  linearsolver.Solve(topology, vector_x, vector_b);

  EXPECT_NEAR(vector_x(0, 0), 0.333333333333333, 1e-06);
  EXPECT_NEAR(vector_x(1, 0), 0.333333333333333, 1e-06);
}

TEST_F(NonlinearSolverTest, primalvariable)
{
  MIRCO::NonLinearSolver nonlinearsolver;
  nonlinearsolver.NonlinearSolve(matrix_, b_vector_, x_vector_, w_, y_);

  EXPECT_NEAR(y_(0, 0), 163213.374921086, 1e-06);
  EXPECT_NEAR(y_(1, 0), 43877.9231473546, 1e-06);
  EXPECT_NEAR(y_(2, 0), 163702.923578063, 1e-06);
  EXPECT_NEAR(y_(3, 0), 55159.5440853170, 1e-06);
  EXPECT_NEAR(y_(4, 0), 10542.1713862417, 1e-06);
  EXPECT_NEAR(y_(5, 0), 53809.0897795325, 1e-06);
  EXPECT_NEAR(y_(6, 0), 148773.412150208, 1e-06);
  EXPECT_NEAR(y_(7, 0), 83711.5732276221, 1e-06);
  EXPECT_NEAR(y_(8, 0), 149262.960807186, 1e-06);
}

TEST_F(NonlinearSolverTest, dualvariable)
{
  MIRCO::NonLinearSolver nonlinearsolver;
  nonlinearsolver.NonlinearSolve(matrix_, b_vector_, x_vector_, w_, y_);

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
  UTILS::ChangeRelativePath(targetfilename, sourcefilename);
  EXPECT_EQ(targetfilename, "../inputfiles/input.dat");
}

TEST(FilesystemUtils, keepabsolutpath)
{
  std::string targetfilename = "/root_dir/home/user/Input/input.dat";
  std::string sourcefilename = "../inputfiles/sourceinput.json";
  UTILS::ChangeRelativePath(targetfilename, sourcefilename);
  EXPECT_EQ(targetfilename, "/root_dir/home/user/Input/input.dat");
}

TEST(readtopology, RMG)
{
  int resolution = 2;
  float Hurst = 0.1;
  bool rand_seed_flag = false;
  int rmg_seed = 95;
  Epetra_SerialDenseMatrix outsurf;
  int N = pow(2, resolution);
  outsurf.Shape(N + 1, N + 1);
  double LateralLength = 1000;

  MIRCO::Rmg surface(resolution, LateralLength, Hurst, rand_seed_flag, rmg_seed);
  surface.GetSurface(outsurf);

  EXPECT_NEAR(outsurf(0, 0), 3790.26122633529, 1e-03);
  EXPECT_NEAR(outsurf(0, 1), 4871.93366647452, 1e-03);
  EXPECT_NEAR(outsurf(0, 2), 11201.8609342306, 1e-03);
  EXPECT_NEAR(outsurf(0, 3), 7003.46380530129, 1e-03);
  EXPECT_NEAR(outsurf(0, 4), 3790.26122633529, 1e-03);
  EXPECT_NEAR(outsurf(1, 0), 11084.2409740010, 1e-03);
  EXPECT_NEAR(outsurf(1, 1), 11886.6633888590, 1e-03);
  EXPECT_NEAR(outsurf(1, 2), 12556.0127153703, 1e-03);
  EXPECT_NEAR(outsurf(1, 3), 5681.76294587438, 1e-03);
  EXPECT_NEAR(outsurf(1, 4), 3648.34675362529, 1e-03);
  EXPECT_NEAR(outsurf(2, 0), 6304.07980664345, 1e-03);
  EXPECT_NEAR(outsurf(2, 1), 3094.97565086270, 1e-03);
  EXPECT_NEAR(outsurf(2, 2), 12745.7212245737, 1e-03);
  EXPECT_NEAR(outsurf(2, 3), 912.374043566341, 1e-03);
  EXPECT_NEAR(outsurf(2, 4), 6659.99316045450, 1e-03);
  EXPECT_NEAR(outsurf(3, 0), 9527.54077507710, 1e-03);
  EXPECT_NEAR(outsurf(3, 1), 3419.42424749242, 1e-03);
  EXPECT_NEAR(outsurf(3, 2), 8800.58346244825, 1e-03);
  EXPECT_NEAR(outsurf(3, 3), 4511.67889999766, 1e-03);
  EXPECT_NEAR(outsurf(3, 4), 1083.02457045383, 1e-03);
  EXPECT_NEAR(outsurf(4, 0), 3790.26122633529, 1e-03);
  EXPECT_NEAR(outsurf(4, 1), 0, 1e-03);
  EXPECT_NEAR(outsurf(4, 2), 4938.79935963518, 1e-03);
  EXPECT_NEAR(outsurf(4, 3), 5669.91464185705, 1e-03);
  EXPECT_NEAR(outsurf(4, 4), 3790.26122633529, 1e-03);
}

TEST(warmstarting, warmstart)
{
  Epetra_SerialDenseMatrix xv0, yv0, xvf, yvf, pf, x0;

  xv0.Shape(1, 3);
  yv0.Shape(1, 3);
  x0.Shape(3, 1);
  xvf.Shape(1, 2);
  yvf.Shape(1, 2);
  pf.Shape(1, 2);

  xv0(0, 0) = 1;
  xv0(0, 1) = 3;
  xv0(0, 2) = 5;

  yv0(0, 0) = 2;
  yv0(0, 1) = 4;
  yv0(0, 2) = 6;

  xvf(0, 0) = 1;
  xvf(0, 1) = 5;

  yvf(0, 0) = 2;
  yvf(0, 1) = 6;

  pf(0, 0) = 10;
  pf(0, 1) = 30;

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
