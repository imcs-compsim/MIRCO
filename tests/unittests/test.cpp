#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialSymDenseMatrix.h>
#include <gtest/gtest.h>
#include <stdlib.h>
#include <vector>
#include "../../src/filesystem_utils.h"
#include "../../src/linearsolver.h"
#include "../../src/nonlinearsolver.h"
#include "nonlinear_solver_test.h"
#include "../../src/topology.h"
#include <memory>

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
  LinearSolver linearsolver;
  linearsolver.Solve(topology, vector_x, vector_b);

  EXPECT_NEAR(vector_x(0, 0), 0.333333333333333, 1e-06);
  EXPECT_NEAR(vector_x(1, 0), 0.333333333333333, 1e-06);
}

TEST_F(NonlinearSolverTest, primalvariable)
{
  NonLinearSolver nonlinearsolver;
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
  NonLinearSolver nonlinearsolver;
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

TEST(readtopology, readfile)
{
  int resolution = 2;
  Epetra_SerialDenseMatrix outsurf;
  std::string filepath = "../bem/Input/sup2.dat";
  ReadFile surface(resolution, filepath);
  surface.GetSurface(outsurf);

  EXPECT_NEAR(outsurf(0, 0), 5.7299175e+01, 1e-06);
  EXPECT_NEAR(outsurf(0, 1), 2.8085118e+01, 1e-06);
  EXPECT_NEAR(outsurf(0, 2), 5.3606249e+01, 1e-06);
  EXPECT_NEAR(outsurf(0, 3), 4.1267492e+01, 1e-06);
  EXPECT_NEAR(outsurf(0, 4), 5.7299175e+01, 1e-06);
  EXPECT_NEAR(outsurf(1, 0), 5.4725631e+01, 1e-06);
  EXPECT_NEAR(outsurf(1, 1), 5.1746180e+01, 1e-06);
  EXPECT_NEAR(outsurf(1, 2), 4.4419276e+01, 1e-06);
  EXPECT_NEAR(outsurf(1, 3), 6.0368219e+01, 1e-06);
  EXPECT_NEAR(outsurf(1, 4), 4.0305024e+01, 1e-06);
  EXPECT_NEAR(outsurf(2, 0), 3.2635027e+01, 1e-06);
  EXPECT_NEAR(outsurf(2, 1), 6.6587744e+01, 1e-06);
  EXPECT_NEAR(outsurf(2, 2), 7.4298269e+01, 1e-06);
  EXPECT_NEAR(outsurf(2, 3), 1.0839591e+02, 1e-06);
  EXPECT_NEAR(outsurf(2, 4), 1.0653673e+02, 1e-06);
  EXPECT_NEAR(outsurf(3, 0), 2.8225470e+01, 1e-06);
  EXPECT_NEAR(outsurf(3, 1), 0.0000000e+00, 1e-06);
  EXPECT_NEAR(outsurf(3, 2), 4.5480444e+01, 1e-06);
  EXPECT_NEAR(outsurf(3, 3), 1.0201105e+02, 1e-06);
  EXPECT_NEAR(outsurf(3, 4), 8.5512730e+01, 1e-06);
  EXPECT_NEAR(outsurf(4, 0), 5.7299175e+01, 1e-06);
  EXPECT_NEAR(outsurf(4, 1), 5.7150667e+01, 1e-06);
  EXPECT_NEAR(outsurf(4, 2), 5.1100808e+01, 1e-06);
  EXPECT_NEAR(outsurf(4, 3), 9.8243100e+01, 1e-06);
  EXPECT_NEAR(outsurf(4, 4), 5.7299175e+01, 1e-06);
}

TEST(readtopology, RMG)
{
  int resolution = 2;
  float Hurst = 0.1;
  bool rand_seed_flag = false;
  Epetra_SerialDenseMatrix outsurf;
  int N = pow(2, resolution);
  outsurf.Shape(N + 1, N + 1);
  
  Rmg surface(resolution, Hurst, rand_seed_flag);
  surface.GetSurface(outsurf);

  EXPECT_NEAR(outsurf(0, 0), 28.2215338276376, 1e-03);
  EXPECT_NEAR(outsurf(0, 1), 36.2753295855912, 1e-03);
  EXPECT_NEAR(outsurf(0, 2), 83.406827899629, 1e-03);
  EXPECT_NEAR(outsurf(0, 3), 52.1463954928658, 1e-03);
  EXPECT_NEAR(outsurf(0, 4), 28.2215338276376, 1e-03);
  EXPECT_NEAR(outsurf(1, 0), 82.5311192966684, 1e-03);
  EXPECT_NEAR(outsurf(1, 1), 88.5057522645388, 1e-03);
  EXPECT_NEAR(outsurf(1, 2), 93.4896404218242, 1e-03);
  EXPECT_NEAR(outsurf(1, 3), 42.3052621232221, 1e-03);
  EXPECT_NEAR(outsurf(1, 4), 27.1648479422294, 1e-03);
  EXPECT_NEAR(outsurf(2, 0), 46.9389540667996, 1e-03);
  EXPECT_NEAR(outsurf(2, 1), 23.0445730666838, 1e-03);
  EXPECT_NEAR(outsurf(2, 2), 94.9021951336941, 1e-03);
  EXPECT_NEAR(outsurf(2, 3), 6.79338672128661, 1e-03);
  EXPECT_NEAR(outsurf(2, 4), 49.5890138641246, 1e-03);
  EXPECT_NEAR(outsurf(3, 0), 70.9402900749234, 1e-03);
  EXPECT_NEAR(outsurf(3, 1), 25.4602528043112, 1e-03);
  EXPECT_NEAR(outsurf(3, 2), 65.527428831836, 1e-03);
  EXPECT_NEAR(outsurf(3, 3), 33.5930588570828, 1e-03);
  EXPECT_NEAR(outsurf(3, 4), 8.064026613556, 1e-03);
  EXPECT_NEAR(outsurf(4, 0), 28.2215338276376, 1e-03);
  EXPECT_NEAR(outsurf(4, 1), 0, 1e-03);
  EXPECT_NEAR(outsurf(4, 2), 36.7733127073012, 1e-03);
  EXPECT_NEAR(outsurf(4, 3), 42.2170752636335, 1e-03);
  EXPECT_NEAR(outsurf(4, 4), 28.2215338276376, 1e-03);
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
