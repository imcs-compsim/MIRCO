#include "contactpredictors.h"
#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>
#include "warmstart.h"

void ContactSetPredictor(int &n0, std::vector<double> &xv0, std::vector<double> &yv0,
    std::vector<double> &b0, double zmax, double Delta, double w_el, std::vector<double> x,
    Epetra_SerialDenseMatrix topology)
{
  std::vector<int> col, row;
  // [ind1,ind2]=find(z>=(zmax-(Delta(s)+w_el(k))));
  double value = zmax - Delta - w_el;
  row.clear();
  col.clear();

  // Data is even, guided makes more sense
  for (int i = 0; i < topology.N(); i++)
  {
    for (int j = 0; j < topology.N(); j++)
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
      xv0[b] = x[col[b]];
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
      yv0[b] = x[row[b]];
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

void InitialGuessPredictor(bool flagwarm, int k, int n0, int nf, std::vector<double> xv0,
    std::vector<double> yv0, std::vector<double> pf, Epetra_SerialDenseMatrix &x0,
    std::vector<double> &b0, std::vector<double> xvf, std::vector<double> yvf)
{
  Epetra_SerialDenseMatrix xv0t, yv0t, xvft, yvft, pft;  // Temporary variables for warmup
  if (flagwarm == 1 && k > 0)
  {
    // x0=warm_x(xv0(1:n0(k),k),yv0(1:n0(k),k),xvf(1:nf(k-1),k-1),yvf(1:nf(k-1),k-1),pf(1:nf(k-1),k-1));
    xv0t.Shape(1, n0);
    yv0t.Shape(1, n0);
    xvft.Shape(1, nf);
    yvft.Shape(1, nf);
    pft.Shape(1, nf);

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

    Warmstart(x0, xv0t, yv0t, xvft, yvft, pft);

    // #pragma omp parallel for schedule(static, 16)  // Always same workload -> Static
    //     for (int i = 0; i < x0temp.M(); i++)
    //     {
    //       x0[i] = x0temp(i, 0);
    //     }
  }
  else
  {
    if (b0.size() > 0)
    {
      x0.Shape(n0, 1);
    }
  }
}