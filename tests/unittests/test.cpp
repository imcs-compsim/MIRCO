#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialSymDenseMatrix.h>
#include <gtest/gtest.h>
#include <stdlib.h>
#include <vector>
#include "../../src/filesystem_utils.h"
#include "../../src/linearsolver.h"
#include "../../src/nonlinearsolver.h"
#include "nonlinear_solver_test.h"

TEST(linearsolver, solves) {
  int systemsize = 2;

  // Build the matrix
  Epetra_SerialSymDenseMatrix topology;
  topology.Shape(systemsize);
  for (int i = 0; i < systemsize; i++) {
    topology(i, i) = 2;
    for (int j = 0; j < i; j++) {
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
  for (int i = 0; i < systemsize; i++) {
    vector_b(i, 0) = 1;
  }

  // Call linear solver
  LinearSolver linearsolver;
  linearsolver.Solve(topology, vector_x, vector_b);

  EXPECT_NEAR(vector_x(0, 0), 0.333333333333333, 1e-06);
  EXPECT_NEAR(vector_x(1, 0), 0.333333333333333, 1e-06);
}

TEST_F(NonlinearSolverTest, primalvariable) {
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

TEST_F(NonlinearSolverTest, dualvariable) {
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

TEST(FilesystemUtils, createrelativepath) {
  std::string targetfilename = "input.dat";
  std::string sourcefilename = "../inputfiles/sourceinput.json";
  ChangeRelativePath(targetfilename, sourcefilename);
  EXPECT_EQ(targetfilename, "../inputfiles/input.dat");
}

TEST(FilesystemUtils, keepabsolutpath) {
  std::string targetfilename = "/root_dir/home/user/Input/input.dat";
  std::string sourcefilename = "../inputfiles/sourceinput.json";
  ChangeRelativePath(targetfilename, sourcefilename);
  EXPECT_EQ(targetfilename, "/root_dir/home/user/Input/input.dat");
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
