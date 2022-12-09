#ifndef SRC_MATRIXSETUP_H_
#define SRC_MATRIXSETUP_H_

#include <Teuchos_SerialDenseMatrix.hpp>
#include <vector>

namespace MIRCO
{
  class MatrixGeneration
  {
   public:
    /**
     * @brief The aim of this function is to create the influence coefficient matrix (Discrete
     * version of Green function)
     *
     * @param A Influence coefficient matrix (Discrete version of Green Function)
     * @param xv0 x-coordinates of the points in contact in the previous iteration.
     * @param yv0 y-coordinates of the points in contact in the previous iteration.
     * @param GridSize Grid size
     * @param CompositeYoungs The composite Young's modulus
     * @param systemsize Number of nodes predicted to be in contact
     */
    void SetUpMatrix(Teuchos::SerialDenseMatrix<int,double>& A, std::vector<double> xv0, std::vector<double> yv0,
        double GridSize, double CompositeYoungs, int systemsize);
    MatrixGeneration() = default;
  };
}  // namespace MIRCO

#endif  // SRC_MATRIXSETUP_H_
