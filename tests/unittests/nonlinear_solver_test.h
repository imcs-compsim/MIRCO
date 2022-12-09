#ifndef NONLINEAR_SOLVER_TEST_H_
#define NONLINEAR_SOLVER_TEST_H_

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <gtest/gtest.h>
#include <vector>

class NonlinearSolverTest : public ::testing::Test
{
 protected:
  NonlinearSolverTest();

  Teuchos::SerialDenseMatrix<int,double> matrix_, b_vector_, w_, y_;
  Teuchos::SerialDenseMatrix<int,double> x_vector_;
};

#endif /* NONLINEAR_SOLVER_TEST_H_ */
