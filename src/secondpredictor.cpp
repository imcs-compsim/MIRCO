#include "secondpredictor.h"
#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>
#include "warmstart.h"

void SecondPredictor(bool flagwarm, int k, int n0, int nf2, std::vector<double> xv0,
    std::vector<double> yv0, std::vector<double> pf, std::vector<double> &x0,
    std::vector<double> &b0, std::vector<double> xvf, std::vector<double> yvf)
{
  Epetra_SerialDenseMatrix xv0t, yv0t, xvft, yvft, pft, x0temp;  // Temporary variables for warmup
  if (flagwarm == 1 && k > 1)
  {
    // x0=warm_x(xv0(1:n0(k),k),yv0(1:n0(k),k),xvf(1:nf(k-1),k-1),yvf(1:nf(k-1),k-1),pf(1:nf(k-1),k-1));
    xv0t.Shape(1, n0);
    yv0t.Shape(1, n0);
    xvft.Shape(1, nf2);
    yvft.Shape(1, nf2);
    pft.Shape(1, nf2);

#pragma omp parallel for schedule(static, 16)  // Always same workload -> Static
    for (int i = 0; i < n0; i++)
    {
      xv0t(0, i) = xv0[i];
      yv0t(0, i) = yv0[i];
    }

#pragma omp parallel for schedule(static, 16)  // Always same workload -> Static
    for (int j = 0; j < nf2; j++)
    {
      xvft(0, j) = xvf[j];
      yvft(0, j) = yvf[j];
      pft(0, j) = pf[j];
    }

    Warmstarter warm1;
    x0temp = warm1.Warmstart2(xv0t, yv0t, xvft, yvft, pft);

#pragma omp parallel for schedule(static, 16)  // Always same workload -> Static
    for (int i = 0; i < x0temp.N(); i++)
    {
      x0[i] = x0temp(i, 1);
    }
  }
  else
  {
    if (b0.size() > 0)
    {
      // x0.Shape(b0.N(), 1);
      x0.resize(b0.size());
    }
  }
}