#ifndef SRC_WARMSTART_H_
#define SRC_WARMSTART_H_

#include <Teuchos_SerialDenseMatrix.hpp>

namespace MIRCO
{
  /**
   * @brief This function is used to determine the nodes which were in contact in the last iteration
   * from the current contact set. This helps in making an initial guess of the nodes in contact in
   * the current iteration and speeds up the computation.
   *
   * @param x0 contact forces at (xvf,yvf) predicted in the previous iteration but are a part of
   * currect predicted contact set.
   * @param xv0 x-coordinates of the points in contact in the previous iteration.
   * @param yv0 y-coordinates of the points in contact in the previous iteration.
   * @param xvf x-coordinates of the points in contact in the previous iteration.
   * @param yvf y-coordinates of the points in contact in the previous iteration.
   * @param pf Contact force at (xvf,yvf) predicted in the previous iteration.
   */
  void Warmstart(Teuchos::SerialDenseMatrix<int,double>& x0, Teuchos::SerialDenseMatrix<int,double> xv0,
      Teuchos::SerialDenseMatrix<int,double> yv0, Teuchos::SerialDenseMatrix<int,double>& xvf, Teuchos::SerialDenseMatrix<int,double>& yvf,
      Teuchos::SerialDenseMatrix<int,double>& pf);
}  // namespace MIRCO

#endif  // SRC_WARMSTART_H_
