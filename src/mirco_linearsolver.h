#ifndef SRC_LINEARSOLVER_H_
#define SRC_LINEARSOLVER_H_

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialSymDenseMatrix.hpp>

namespace MIRCO
{
  namespace LinearSolver
  {
    /**
     * @brief Solve the linear system matrix * vector_x = vector_b
     *
     * @param matrix Influence coefficient matrix (Discrete version of Green Function)
     * @param vector_b RHS
     *
     * @return vector_x Solution to the linear system
     */
    Teuchos::SerialDenseVector<int, double> Solve(
        Teuchos::SerialSymDenseMatrix<int, double>& matrix,
        Teuchos::SerialDenseVector<int, double>& vector_b);
  };  // namespace LinearSolver
}  // namespace MIRCO

#endif  // SRC_LINEARSOLVER_H_
