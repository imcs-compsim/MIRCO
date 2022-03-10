#ifndef SRC_WARMSTART_H_
#define SRC_WARMSTART_H_

#include <Epetra_SerialSymDenseMatrix.h>

void Warmstart(Epetra_SerialDenseMatrix& x0, Epetra_SerialDenseMatrix xv0,
    Epetra_SerialDenseMatrix yv0, Epetra_SerialDenseMatrix& xvf, Epetra_SerialDenseMatrix& yvf,
    Epetra_SerialDenseMatrix& pf);

#endif  // SRC_WARMSTART_H_
