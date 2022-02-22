#include "contactsetpredictor.h"
#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>

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