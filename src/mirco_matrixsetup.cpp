#include "mirco_matrixsetup.h"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <vector>

void MIRCO::MatrixGeneration::SetUpMatrix(Teuchos::SerialDenseMatrix<int, double>& A,
    std::vector<double> xv0, std::vector<double> yv0, double GridSize, double CompositeYoungs,
    double CompositePoissonsRatio, int systemsize, bool PressureGreenFunFlag)
{
  double pi = atan(1) * 4;
  if (PressureGreenFunFlag)
  {
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
    A.scale((1 - pow(CompositePoissonsRatio, 2)) / (pi * CompositeYoungs));
  }
  else
  {
    double r, raggio = GridSize / 2;
    double C = 1 / (CompositeYoungs * pi * raggio);

#pragma omp parallel for schedule(static, 16)  // Always same workload -> static
    for (int i = 0; i < systemsize; i++)
    {
      A(i, i) = 1 * C;
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
