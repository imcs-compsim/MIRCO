#ifndef SRC_LINEARSOLVER_H_
#define SRC_LINEARSOLVER_H_

#include <Epetra_SerialSpdDenseSolver.h>
#include <Epetra_SerialSymDenseMatrix.h>

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
    void Solve(Epetra_SerialSymDenseMatrix& matrix, Epetra_SerialDenseMatrix& vector_x,
        Epetra_SerialDenseMatrix& vector_b);
    LinearSolver() = default;
  };
}  // namespace MIRCO

#endif  // SRC_LINEARSOLVER_H_
