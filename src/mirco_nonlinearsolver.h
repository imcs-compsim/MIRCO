#ifndef SRC_NONLINEARSOLVER_H_
#define SRC_NONLINEARSOLVER_H_

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <vector>

namespace MIRCO
{
  namespace NonLinearSolver
  {
    /**
     * @brief Solve the non-linear set of equations using Non-Negative Least Squares (NNLS) method
     * and return the contact force solution y.
     *
     * @param matrix Influence coefficient matrix (Discrete version of Green Function)
     * @param b0 Indentation value of the half space at the predicted points of contact.
     * @param y0 contact forces at (xvf,yvf) predicted in the previous iteration but are a part of
     * currect predicted contact set.
     * @param w Gap between the point on the topology and the half space
     *
     * @return y Solution containing force
     */
    Teuchos::SerialDenseVector<int, double> Solve(
        const Teuchos::SerialDenseMatrix<int, double>& matrix, const std::vector<double>& b0,
        const Teuchos::SerialDenseMatrix<int, double>& y0,
        Teuchos::SerialDenseMatrix<int, double>& w);
  };  // namespace NonLinearSolver
}  // namespace MIRCO

#endif  // SRC_NONLINEARSOLVER_H_
