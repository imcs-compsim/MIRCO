#ifndef SRC_COMPUTERESIDUAL_H_
#define SRC_COMPUTERESIDUAL_H_

#include <Epetra_SerialSymDenseMatrix.h>

void ComputeResidual(Epetra_SerialDenseMatrix A, Epetra_SerialDenseMatrix y,
    Epetra_SerialDenseMatrix b0new, Epetra_SerialDenseMatrix w, Epetra_SerialDenseMatrix &res1);

#endif  // SRC_COMPUTERESIDUAL_H_