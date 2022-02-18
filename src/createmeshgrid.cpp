#include "createmeshgrid.h"
#include <cmath>
#include <vector>

void CreateMeshgrid(std::vector<double> &x, int iter, double delta)
{
#pragma omp parallel for schedule(static, 16)  // Same amount of work -> static
  for (int i = 0; i < iter; i++)
  {
    x[i] = (delta / 2) + i * delta;
  }
}