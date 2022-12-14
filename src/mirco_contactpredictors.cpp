#include "mirco_contactpredictors.h"
#include <Teuchos_SerialDenseMatrix.hpp>
#include <vector>
#include "mirco_warmstart.h"

void MIRCO::ContactSetPredictor(int &n0, std::vector<double> &xv0, std::vector<double> &yv0,
    std::vector<double> &b0, double zmax, double Delta, double w_el, std::vector<double> &meshgrid,
    Teuchos::SerialDenseMatrix<int, double> &topology)
{
  std::vector<int> col, row;
  double value = zmax - Delta - w_el;
  row.clear();
  col.clear();

  // Data is even, guided makes more sense
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

#pragma omp for schedule(guided, 16)  // Always same workload but testing might be good -> Guided?
  for (int b = 0; b < n0; b++)
  {
    try
    {
      xv0[b] = meshgrid[col[b]];
    }
    catch (const std::exception &e)
    {
    }
  }

#pragma omp parallel for schedule(guided, 16)  // Same
  for (int b = 0; b < n0; b++)
  {
    try
    {
      yv0[b] = meshgrid[row[b]];
    }
    catch (const std::exception &e)
    {
    }
  }

#pragma omp parallel for schedule(guided, 16)  // Same
  for (int b = 0; b < n0; b++)
  {
    try
    {
      b0[b] = Delta + w_el - (zmax - topology(row[b], col[b]));
    }
    catch (const std::exception &e)
    {
    }
  }
}

void MIRCO::InitialGuessPredictor(bool WarmStartingFlag, int k, int n0, int nf,
    std::vector<double> xv0, std::vector<double> yv0, std::vector<double> pf,
    Teuchos::SerialDenseMatrix<int, double> &x0, std::vector<double> &b0, std::vector<double> xvf,
    std::vector<double> yvf)
{
  Teuchos::SerialDenseMatrix<int, double> xv0t, yv0t, xvft, yvft,
      pft;  // Temporary variables for warmup
  if (WarmStartingFlag == 1 && k > 0)
  {
    xv0t.shape(1, n0);
    yv0t.shape(1, n0);
    xvft.shape(1, nf);
    yvft.shape(1, nf);
    pft.shape(1, nf);

#pragma omp parallel for schedule(static, 16)  // Always same workload -> Static
    for (int i = 0; i < n0; i++)
    {
      xv0t(0, i) = xv0[i];
      yv0t(0, i) = yv0[i];
    }

#pragma omp parallel for schedule(static, 16)  // Always same workload -> Static
    for (int j = 0; j < nf; j++)
    {
      xvft(0, j) = xvf[j];
      yvft(0, j) = yvf[j];
      pft(0, j) = pf[j];
    }

    MIRCO::Warmstart(x0, xv0t, yv0t, xvft, yvft, pft);
  }
  else
  {
    if (b0.size() > 0)
    {
      x0.shape(n0, 1);
    }
  }
}
