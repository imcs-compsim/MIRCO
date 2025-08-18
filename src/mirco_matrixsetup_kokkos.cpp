#include "mirco_matrixsetup_kokkos.h"

#include <math.h>

namespace MIRCO
{
  ViewMatrix_d MatrixGeneration::SetupMatrix(const ViewVector_d xv0, const ViewVector_d yv0,
      const double GridSize, const double CompositeYoungs, const int systemsize,
      const bool PressureGreenFunFlag)
  {
    constexpr double pi = M_PI;

    const double frac_GridSize_2 = GridSize / 2;

    ViewMatrix_d H_d("MatrixGeneration::SetupMatrix(); H_d", xv0.size(), xv0.size());
    if (PressureGreenFunFlag)
    {
      // The pressure-based Green's function is based on the work of Pohrt and Li (2014)
      // https://doi.org/10.1134/S1029959914040109
      // Please look at equation 12 of the paper mentioned above.
      // ((1-nu)/2*pi*G) from the equation is replaced with (1/pi*CompositeYoungs) here.
      // The paper uses a decoupled shear modulus and Poisson's ratio. We use a composite Young's
      // modulus here, instead.

      // Note: KOKKOS_LAMBDA will automatically capture const variables from the outer scope into
      // device-space from host-space (but not non-const variables!)
      Kokkos::parallel_for(
          Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {systemsize, systemsize}),
          KOKKOS_LAMBDA(const int i, const int j) {
            double k = xv0(i) - xv0(j) + frac_GridSize_2;
            double l = xv0(i) - xv0(j) - frac_GridSize_2;
            double m = yv0(i) - yv0(j) + frac_GridSize_2;
            double n = yv0(i) - yv0(j) - frac_GridSize_2;

            H_d(i, j) = 1 / (pi * CompositeYoungs) *
                        (k * log((sqrt(k * k + m * m) + m) / (sqrt(k * k + n * n) + n)) +
                            l * log((sqrt(l * l + n * n) + n) / (sqrt(l * l + m * m) + m)) +
                            m * log((sqrt(m * m + k * k) + k) / (sqrt(m * m + l * l) + l)) +
                            n * log((sqrt(n * n + l * l) + l) / (sqrt(n * n + k * k) + k)));
          });
    }

    else
    {
      const double C = 1 / (CompositeYoungs * pi * frac_GridSize_2);

      // TODO: For better effiency, try using teams instead of MDRangePolicy
      Kokkos::parallel_for(
          Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {systemsize, systemsize}),
          KOKKOS_LAMBDA(const int i, const int j) {
            if (j >= i) return;
            const double tmp1 = xv0(j) - xv0(i);
            const double tmp2 = yv0(j) - yv0(i);
            const double r = sqrt(tmp1 * tmp1 + tmp2 * tmp2);
            const double tmp3 = C * asin(frac_GridSize_2 / r);
            H_d(i, j) = tmp3;
            H_d(j, i) = tmp3;
          });

      Kokkos::parallel_for(systemsize, KOKKOS_LAMBDA(const int i) { H_d(i, i) = C; });
    }

    return H_d;
  }

}  // namespace MIRCO
