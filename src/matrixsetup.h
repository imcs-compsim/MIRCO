#ifndef SRC_MATRIXSETUP_H_
#define SRC_MATRIXSETUP_H_

#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>

namespace MIRCO
{
  class MatrixGeneration
  {
   public:
    void SetUpMatrix(Epetra_SerialDenseMatrix& A, std::vector<double> xv0, std::vector<double> yv0,
        double delta, double E, int systemsize);
    MatrixGeneration() = default;
  };
}  // namespace MIRCO

#endif  // SRC_MATRIXSETUP_H_
