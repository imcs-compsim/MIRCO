#ifndef NONLINEAR_SOLVER_TEST_H_
#define NONLINEAR_SOLVER_TEST_H_

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialSymDenseMatrix.h>
#include <gtest/gtest.h>
#include <vector>

class NonlinearSolverTest : public ::testing::Test
{
 protected:
  NonlinearSolverTest();

  Epetra_SerialDenseMatrix matrix_, b_vector_, w_, y_;
  Epetra_SerialDenseMatrix x_vector_;
};

#endif /* NONLINEAR_SOLVER_TEST_H_ */
