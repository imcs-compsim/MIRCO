#include "computeresidual.h"
#include <Epetra_SerialSymDenseMatrix.h>

void ComputeResidual(Epetra_SerialDenseMatrix A, Epetra_SerialDenseMatrix y,
    Epetra_SerialDenseMatrix b0new, Epetra_SerialDenseMatrix w, Epetra_SerialDenseMatrix &res1)
{
  if (A.M() != y.N())
  {
    std::runtime_error("Error 1: Matrix dimensions imcompatible");
  }
  res1.Shape(A.M(), A.M());

  // res1=A*sol-b0(:,k)-wsol;
#pragma omp parallel for schedule(static, 16)  // Always same workload -> Static
  for (int x = 0; x < A.N(); x++)
  {
    for (int z = 0; z < y.M(); z++)
    {
      res1(x, 0) += A(x, z) * y(z, 0);
    }
    res1(x, 0) -= b0new(x, 0) + w(x, 0);  // [...]-b0(:,k) - wsol;
  }
}