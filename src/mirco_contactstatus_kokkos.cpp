#include "mirco_contactstatus_kokkos.h"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <cmath>
#include <vector>

void MIRCO::ComputeContactNodes(std::vector<double> &xvf, std::vector<double> &yvf,
    std::vector<double> &pf, int &nf, const Teuchos::SerialDenseMatrix<int, double> y,
    const std::vector<double> xv0, const std::vector<double> yv0)
{
  xvf.clear();
  xvf.resize(y.numRows());
  yvf.clear();
  yvf.resize(y.numRows());
  pf.resize(y.numRows());
  nf = 0;
  // @} Parallelizing this slows down program, so removed it.

#pragma omp for schedule(guided, 16)
  for (int i = 0; i < y.numRows(); i++)
  {
    if (y(i, 0) != 0)
    {
#pragma omp critical
      {
        xvf[nf] = xv0[i];
        yvf[nf] = yv0[i];
        pf[nf] = y(i, 0);
        nf += 1;
      }
    }
  }
}

void MIRCO::ComputeContactForceAndArea(std::vector<double> &force0, std::vector<double> &area0, int nf, std::vector<double> pf, int k, double GridSize, double LateralLength, bool PressureGreenFunFlag)
{
  force0.push_back(0);
  double sum = 0;
#pragma omp parallel for schedule(static, 16) reduction(+ : sum)  // Always same workload -> Static!
  for (int i = 0; i < nf; i++)
  {
    if (PressureGreenFunFlag)
    {
      sum += pf[i] * pow(GridSize, 2);
    }
    else
    {
      sum += pf[i];
    }
  }
  force0[k] += sum;
  area0.push_back(nf * (pow(GridSize, 2) / pow(LateralLength, 2)));
}
