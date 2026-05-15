#ifndef SRC_MATRIXSETUP_H_
#define SRC_MATRIXSETUP_H_

#include "mirco_kokkostypes.h"

namespace MIRCO
{
  /**
   * @brief Create the influence coefficient matrix (Discrete version of Green's function)
   *
   * @param[in] xv0 x-coordinates of the points in contact in the previous iteration.
   * @param[in] yv0 y-coordinates of the points in contact in the previous iteration.
   * @param[in] GridSize Grid size (length of each cell)
   * @param[in] CompositeYoungs The composite Young's modulus
   * @param[in] systemsize Number of nodes predicted to be in contact
   * @param[in] PressureGreenFunFlag Flag to use Green function based on uniform pressure instead
   * of point force
   *
   * @return Influence coefficient matrix (Discrete version of Green Function) (usually denoted H)
   */
  ViewMatrix_d SetupMatrix(const ViewVector_d xv0, const ViewVector_d yv0, const double GridSize,
      const double CompositeYoungs, const int systemsize, const bool PressureGreenFunFlag);

  /**
   * @brief Compute one entry of the full influence coefficient matriix. Use when memory is too
   * constrained to store the full matrix.
   *
   * @param[in] ix x index of first point
   * @param[in] iy y index of first point
   * @param[in] jx x index of second point
   * @param[in] jy y index of second point
   * @param[in] GridSize Grid size (length of each cell)
   * @param[in] CompositeYoungs The composite Young's modulus
   * @param[in] N Element count along one direction
   * @param[in] PressureGreenFunFlag Flag to use Green function based on uniform pressure instead
   * of point force
   *
   * @return Matrix entry of H
   */

  KOKKOS_INLINE_FUNCTION
  double SetupMatrixOneEntry(const int ix, const int iy, const int jx, const int jy,
      const double GridSize, const double CompositeYoungs, const int N,
      const bool PressureGreenFunFlag)
  {
    constexpr double pi = M_PI;
    const double frac_GridSize_2 = GridSize / 2;

    const double xi = ix * GridSize;
    const double yi = iy * GridSize;
    const double xj = jx * GridSize;
    const double yj = jy * GridSize;

    if (PressureGreenFunFlag)
    {
      const double k = xi - xj + frac_GridSize_2;
      const double l = k - GridSize;
      const double m = yi - yj + frac_GridSize_2;
      const double n = m - GridSize;

      return 1.0 / (pi * CompositeYoungs) *
             (k * log((sqrt(k * k + m * m) + m) / (sqrt(k * k + n * n) + n)) +
                 l * log((sqrt(l * l + n * n) + n) / (sqrt(l * l + m * m) + m)) +
                 m * log((sqrt(m * m + k * k) + k) / (sqrt(m * m + l * l) + l)) +
                 n * log((sqrt(n * n + l * l) + l) / (sqrt(n * n + k * k) + k)));
    }

    else
    {
      const double C = 1 / (CompositeYoungs * pi * frac_GridSize_2);

      if (ix == jx && iy == jy)
        return C;

      else
      {
        const double tmp1 = xj - xi;
        const double tmp2 = yj - yi;
        const double r = sqrt(tmp1 * tmp1 + tmp2 * tmp2);
        return C * asin(frac_GridSize_2 / r);
      }
    }
  }
}  // namespace MIRCO

#endif  // SRC_MATRIXSETUP_H_
