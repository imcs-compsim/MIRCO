#include "mirco_warmstart.h"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <vector>

void MIRCO::Warmstart(Teuchos::SerialDenseMatrix<int, double>& x0, const std::vector<double>& xv0,
    const std::vector<double>& yv0, const std::vector<double>& xvf, const std::vector<double>& yvf,
    const std::vector<double>& pf)
{
  x0.shape(xv0.size(), 1);

  for (size_t i = 0; i < xv0.size(); i++)
  {
    auto it_x = std::find(xvf.begin(), xvf.end(), xv0[i]);
    auto it_y = std::find(yvf.begin(), yvf.end(), yv0[i]);

    if (it_x != xvf.end() && it_y != yvf.end() &&
        std::distance(xvf.begin(), it_x) == std::distance(yvf.begin(), it_y))
    {
      size_t j = std::distance(xvf.begin(), it_x);
      x0(i, 0) = pf[j];
    }
  }
}
