#include <Teuchos_SerialDenseMatrix.hpp>
#include <vector>

#include "mirco_matrixsetup.h"

void MIRCO::MatrixGeneration::SetUpMatrix(Teuchos::SerialDenseMatrix<int, double>& A,
    std::vector<double> xv0, std::vector<double> yv0, double GridSize, double CompositeYoungs,
    int systemsize)
{
  double r, pi = atan(1) * 4, raggio = GridSize / 2;
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
