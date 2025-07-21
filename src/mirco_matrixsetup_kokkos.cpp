#include "mirco_matrixsetup_kokkos.h"

#include <math.h>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <vector>

ViewMatrix_d MIRCO::MatrixGeneration::SetupMatrix(const ViewVector_d& xv0, const ViewVector_d& yv0,
    const double GridSize, const double CompositeYoungs, const int systemsize,
    const bool PressureGreenFunFlag)
{
  constexpr double pi = M_PI;

  ViewMatrix_d H_d("H", xv0.size(), xv0.size());
  if (PressureGreenFunFlag)
  {
    // The pressure-based Green's function is based on the work of Pohrt and Li (2014)
    // https://doi.org/10.1134/S1029959914040109
    // Please look at equation 12 of the paper mentioned above.
    // ((1-nu)/2*pi*G) from the equation is replaced with (1/pi*CompositeYoungs) here.
    // The paper uses a decoupled shear modulus and Poisson's ratio. We use a composite Young's
    // modulus here, instead.

    Kokkos::parallel_for(
        "H full", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {systemsize, systemsize}),
        KOKKOS_LAMBDA(const int i, const int j) {
          double k = xv0(i) - xv0(j) + GridSize / 2;
          double l = xv0(i) - xv0(j) - GridSize / 2;
          double m = yv0(i) - yv0(j) + GridSize / 2;
          double n = yv0(i) - yv0(j) - GridSize / 2;

          H_d(i, j) = 1 / (pi * CompositeYoungs) *
                      (k * log((sqrt(k * k + m * m) + m) / (sqrt(k * k + n * n) + n)) +
                          l * log((sqrt(l * l + n * n) + n) / (sqrt(l * l + m * m) + m)) +
                          m * log((sqrt(m * m + k * k) + k) / (sqrt(m * m + l * l) + l)) +
                          n * log((sqrt(n * n + l * l) + l) / (sqrt(n * n + k * k) + k)));
        });
  }

  else
  {
    double raggio = GridSize / 2;
    double C = 1 / (CompositeYoungs * pi * raggio);
    Kokkos::parallel_for(
        "H diag", systemsize,
        KOKKOS_LAMBDA(
            const int i) {  // Kokkos::RangePolicy(0, systemsize), KOKKOS_LAMBDA(const int i) {
          H_d(i, i) = C;
        });

    Kokkos::parallel_for(
        "H off-diag", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {systemsize, systemsize}),
        KOKKOS_LAMBDA(const int i, const int j) {
          double r = sqrt(pow((xv0(j) - xv0(i)), 2) + pow((yv0(j) - yv0(i)), 2));
          H_d(i, j) = C * asin(raggio / r);
          H_d(j, i) = C * asin(raggio / r);
        });
  }

  return H_d;
}
