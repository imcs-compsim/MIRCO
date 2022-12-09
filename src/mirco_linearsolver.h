#ifndef SRC_LINEARSOLVER_H_
#define SRC_LINEARSOLVER_H_

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialSymDenseMatrix.hpp>

namespace MIRCO
{
  class LinearSolver
  {
   public:
    /**
     * @brief A linear solver to solve matrix*vector_x=vector_b
     *
     * @param matrix Influence coefficient matrix (Discrete version of Green Function)
     * @param vector_x Solution
     * @param vector_b RHS
     */
    void Solve(Teuchos::SerialSymDenseMatrix<int,double>& matrix, Teuchos::SerialDenseMatrix<int,double>& vector_x,
        Teuchos::SerialDenseMatrix<int,double>& vector_b);
    LinearSolver() = default;
  };
}  // namespace MIRCO

#endif  // SRC_LINEARSOLVER_H_
