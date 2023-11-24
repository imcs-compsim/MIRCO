#ifndef NONLINEAR_SOLVER_TEST_H_
#define NONLINEAR_SOLVER_TEST_H_

#include <gtest/gtest.h>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <vector>

class NonlinearSolverTest : public ::testing::Test
{
 protected:
  NonlinearSolverTest();

  Teuchos::SerialDenseMatrix<int, double> matrix_, w_;
  Teuchos::SerialDenseVector<int, double> x_vector_, y_;
  std::vector<double> b_vector_;
};

#endif /* NONLINEAR_SOLVER_TEST_H_ */
