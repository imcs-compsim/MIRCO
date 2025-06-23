#include "mirco_matrixsetup.h"

#include <math.h>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <vector>

void MIRCO::MatrixGeneration::SetUpMatrix(Teuchos::SerialDenseMatrix<int, double>& A,
    const std::vector<double>& xv0, const std::vector<double>& yv0, const double GridSize,
    const double CompositeYoungs, const int systemsize, const bool PressureGreenFunFlag)
{
  double pi = M_PI;
  A.shape(xv0.size(), xv0.size());
  if (PressureGreenFunFlag)
  {
    // The pressure-based Green's function is based on the work of Pohrt and Li (2014)
    // https://doi.org/10.1134/S1029959914040109
    // Please look at equation 12 of the paper mentioned above.
    // ((1-nu)/2*pi*G) from the equation is replaced with (1/pi*CompositeYoungs) here.
    // The paper uses a decoupled shear modulus and Poisson's ratio. We use a composite Young's
    // modulus here, instead.
    for (int i = 0; i < systemsize; i++)
    {
      for (int j = 0; j < systemsize; j++)
      {
        double k = xv0[i] - xv0[j] + GridSize / 2;
        double l = xv0[i] - xv0[j] - GridSize / 2;
        double m = yv0[i] - yv0[j] + GridSize / 2;
        double n = yv0[i] - yv0[j] - GridSize / 2;

        A(i, j) = (k * log((sqrt(k * k + m * m) + m) / (sqrt(k * k + n * n) + n)) +
                   l * log((sqrt(l * l + n * n) + n) / (sqrt(l * l + m * m) + m)) +
                   m * log((sqrt(m * m + k * k) + k) / (sqrt(m * m + l * l) + l)) +
                   n * log((sqrt(n * n + l * l) + l) / (sqrt(n * n + k * k) + k)));
      }
    }
    A.scale(1 / (pi * CompositeYoungs));
  }
  else
  {
    double r, raggio = GridSize / 2;
    double C = 1 / (CompositeYoungs * pi * raggio);

#pragma omp parallel for schedule(static, 16)  // Always same workload -> static
    for (int i = 0; i < systemsize; i++)
    {
      A(i, i) = C;
    }

#pragma omp parallel for schedule(static, 16) private(r)  // Always same workload -> static
    // Every iteration needs to have a different r! -> private(r)
    for (int i = 0; i < systemsize; i++)
    {
      for (int j = 0; j < i; j++)
      {
        r = sqrt(pow((xv0[j] - xv0[i]), 2) + pow((yv0[j] - yv0[i]), 2));
        A(i, j) = C * asin(raggio / r);
        A(j, i) = C * asin(raggio / r);
      }
    }
  }
}
