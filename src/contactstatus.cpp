#include "contactstatus.h"
#include <Epetra_SerialSymDenseMatrix.h>
#include <cmath>
#include <vector>

void ComputeContactNodes(std::vector<double> &xvf, std::vector<double> &yvf,
    std::vector<double> &pf, int &cont, int &nf, Epetra_SerialDenseMatrix y,
    std::vector<double> xv0, std::vector<double> yv0)
{
  xvf.clear();
  xvf.resize(y.M());
  yvf.clear();
  yvf.resize(y.M());
  pf.resize(y.M());
  cont = 0;
  // @} Parallelizing this slows down program, so removed it.

#pragma omp for schedule(guided, 16)
  for (int i = 0; i < y.M(); i++)
  {
    if (y(i, 0) != 0)
    {
#pragma omp critical
      {
        xvf[cont] = xv0[i];
        yvf[cont] = yv0[i];
        pf[cont] = y(i, 0);
        cont += 1;
      }
    }
  }

  nf = cont;
}

void ComputeContactForceAndArea(std::vector<double> &force0, std::vector<double> &area0,
    double &w_el, int nf, std::vector<double> pf, int k, double delta, double lato, double k_el)
{
  force0.push_back(0);
  double sum = 0;
#pragma omp parallel for schedule(static, 16) reduction(+ : sum)  // Always same workload -> Static!
  for (int i = 0; i < nf; i++)
  {
    sum += pf[i];
  }
  force0[k] += sum;
  area0.push_back(nf * (pow(delta, 2) / pow(lato, 2)) * 100);
  w_el = force0[k] / k_el;
}