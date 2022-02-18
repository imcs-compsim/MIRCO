#include "computecontactnodes.h"
#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>


void ComputeContactNodes(std::vector<double> &xvf, std::vector<double> &yvf,
    std::vector<double> &pf, int &cont, double &nf, Epetra_SerialDenseMatrix y,
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