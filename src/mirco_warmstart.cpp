#include "mirco_warmstart.h"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <vector>

void MIRCO::Warmstart(Teuchos::SerialDenseMatrix<int, double>& x0,
    Teuchos::SerialDenseMatrix<int, double> xv0, Teuchos::SerialDenseMatrix<int, double> yv0,
    Teuchos::SerialDenseMatrix<int, double>& xvf, Teuchos::SerialDenseMatrix<int, double>& yvf,
    Teuchos::SerialDenseMatrix<int, double>& pf)
{
  x0.shape(xv0.numCols(), 1);

  for (int i = 0; i < xv0.numCols(); i++)
  {
    for (int j = 0; j < xvf.numCols(); j++)
    {
      if (xvf(0, j) == xv0(0, i) && yvf(0, j) == yv0(0, i))
      {
        x0(i, 0) = pf(0, j);
      }
    }
  }
}
