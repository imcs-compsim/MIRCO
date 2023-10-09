#include "mirco_warmstart.h"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <vector>

void MIRCO::Warmstart(Teuchos::SerialDenseMatrix<int, double>& x0, const std::vector<double>& xv0,
    const std::vector<double>& yv0, const std::vector<double>& xvf, const std::vector<double>& yvf,
    const std::vector<double>& pf)
{
  x0.shape(xv0.size(), 1);

  for (long unsigned int i = 0; i < xv0.size(); i++)
  {
    for (long unsigned int j = 0; j < xvf.size(); j++)
    {
      if (xvf[j] == xv0[i] && yvf[j] == yv0[i])
      {
        x0(i, 0) = pf[j];
      }
    }
  }
}
