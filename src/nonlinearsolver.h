#ifndef SRC_NONLINEARSOLVER_H_
#define SRC_NONLINEARSOLVER_H_

#include <Epetra_SerialSpdDenseSolver.h>
#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>

namespace MIRCO
{
  class NonLinearSolver
  {
   public:
    /**
     * @brief /The aim of this function is to solve the non-linear set of equations using
     * Non-Negative Least Squares (NNLS) method.
     *
     * @param matrix Influence coefficient matrix (Discrete version of Green Function)
     * @param b0 Indentation value of the half space at the point of contact.
     * @param y0 contact forces at (xvf,yvf) predicted in the previous iteration but are a part of
     * currect predicted contact set.
     * @param w Defined as (u - u(bar))
     * @param y Solution containing force
     */
    void NonlinearSolve(Epetra_SerialDenseMatrix& matrix, Epetra_SerialDenseMatrix& b0,
        Epetra_SerialDenseMatrix& y0, Epetra_SerialDenseMatrix& w, Epetra_SerialDenseMatrix& y);
    NonLinearSolver() = default;
  };
}  // namespace MIRCO

#endif  // SRC_NONLINEARSOLVER_H_
