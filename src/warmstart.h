#ifndef SRC_WARMSTART_H_
#define SRC_WARMSTART_H_

#include <Epetra_SerialSymDenseMatrix.h>

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
  void Warmstart(Epetra_SerialDenseMatrix& x0, Epetra_SerialDenseMatrix xv0,
      Epetra_SerialDenseMatrix yv0, Epetra_SerialDenseMatrix& xvf, Epetra_SerialDenseMatrix& yvf,
      Epetra_SerialDenseMatrix& pf);
}  // namespace MIRCO

#endif  // SRC_WARMSTART_H_
