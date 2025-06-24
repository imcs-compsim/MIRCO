#ifndef SRC_MATRIXSETUP_H_
#define SRC_MATRIXSETUP_H_

#include <Teuchos_SerialDenseMatrix.hpp>
#include <vector>

namespace MIRCO
{
  namespace MatrixGeneration
  {
    /**
     * @brief The aim of this function is to create the influence coefficient matrix (Discrete
     * version of Green function)
     *
     * @param A Influence coefficient matrix (Discrete version of Green Function)
     * @param xv0 x-coordinates of the points in contact in the previous iteration.
     * @param yv0 y-coordinates of the points in contact in the previous iteration.
     * @param GridSize Grid size (length of each cell)
     * @param CompositeYoungs The composite Young's modulus
     * @param systemsize Number of nodes predicted to be in contact
     * @param PressureGreenFunFlag Flag to use Green function based on uniform pressure instead of
     * point force
     */
    void SetUpMatrix(Teuchos::SerialDenseMatrix<int, double>& A, const std::vector<double>& xv0,
        const std::vector<double>& yv0, const double GridSize, const double CompositeYoungs,
        const int systemsize, const bool PressureGreenFunFlag);
  };  // namespace MatrixGeneration
}  // namespace MIRCO

#endif  // SRC_MATRIXSETUP_H_
