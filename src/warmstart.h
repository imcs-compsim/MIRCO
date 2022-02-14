#ifndef SRC_WARMSTART_H_
#define SRC_WARMSTART_H_

#include <Epetra_SerialSymDenseMatrix.h>

class Warmstarter
{
 public:
  Epetra_SerialDenseMatrix Warmstart2(Epetra_SerialDenseMatrix xv0, Epetra_SerialDenseMatrix yv0,
      Epetra_SerialDenseMatrix &xvf, Epetra_SerialDenseMatrix &yvf, Epetra_SerialDenseMatrix &pf);
  Warmstarter() = default;
};

#endif  // SRC_WARMSTART_H_
