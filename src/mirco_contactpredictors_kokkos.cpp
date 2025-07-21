#include "mirco_contactpredictors_kokkos.h"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <vector>

#include "mirco_warmstart.h"

void MIRCO::ContactSetPredictor(int &n0, std::vector<double> &xv0, std::vector<double> &yv0,
    std::vector<double> &b0, double zmax, double Delta, double w_el,
    const std::vector<double> &meshgrid, const Teuchos::SerialDenseMatrix<int, double> &topology)
{
  std::vector<int> col, row;
  double value = zmax - Delta - w_el;
  // row.clear();//#
  // col.clear();//# completely unecessary

  // Data is even, guided makes more sense//#we parallize
  for (int i = 0; i < topology.numCols(); i++)
  {
    for (int j = 0; j < topology.numCols(); j++)
    {
      if (topology(i, j) >= value)
      {
        row.push_back(i);
        col.push_back(j);
      }
    }
  }

  n0 = col.size();

  // @{
  xv0.clear();
  xv0.resize(n0);
  yv0.clear();
  yv0.resize(n0);
  b0.clear();
  b0.resize(n0);
  // @} Parallelizing slows down program here, so not parallel

#pragma omp parallel for schedule( \
        guided, 16)  // Always same workload but testing might be good -> Guided?
  for (int i = 0; i < n0; i++)
  {
    try
    {
      xv0[i] = meshgrid[col[i]];
      yv0[i] = meshgrid[row[i]];
      b0[i] = Delta + w_el - (zmax - topology(row[i], col[i]));
    }
    catch (const std::exception &e)
    {
    }
  }
}

void MIRCO::InitialGuessPredictor(bool WarmStartingFlag, int k, int n0,
    const std::vector<double> &xv0, const std::vector<double> &yv0, const std::vector<double> &pf,
    Teuchos::SerialDenseMatrix<int, double> &x0, const std::vector<double> &b0,
    const std::vector<double> &xvf, const std::vector<double> &yvf)
{
  if (WarmStartingFlag && k > 0)
  {
    MIRCO::Warmstart(x0, xv0, yv0, xvf, yvf, pf);
  }
  else
  {
    if (b0.size() > 0)
    {
      x0.shape(n0, 1);
    }
  }
}
