#ifndef SRC_MATRIXSETUP_H_
#define SRC_MATRIXSETUP_H_

#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>

class MatrixGeneration
{
 public:
  void SetUpMatrix(Epetra_SerialDenseMatrix& A, std::vector<double> xv0, std::vector<double> yv0,
      double delta, double E, int systemsize);
  MatrixGeneration() = default;
};

#endif  // SRC_MATRIXSETUP_H_
