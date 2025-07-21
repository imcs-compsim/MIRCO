#include "mirco_warmstart_kokkos.h"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <algorithm>
#include <vector>


// # TODO: we now have two overloads, one for host and one for device. compare th performance of
// them including the deep_copy needed for host overload variant.

ViewVector_h MIRCO::Warmstart(const ViewVector_h& xv0, const ViewVector_h& yv0,
    const ViewVector_h& xvf, const ViewVector_h& yvf, const ViewVector_h& pf)
{
  ViewVector_h p0;

  for (size_t i = 0; i < xv0.extent(0); ++i)
  {
    auto it_x = std::find(xvf.data(), xvf.data() + xvf.extent(0), xv0(i));
    auto it_y = std::find(yvf.data(), yvf.data() + yvf.extent(0), yv0(i));

    if (it_x != xvf.data() + xvf.extent(0) && it_y != yvf.data() + yvf.extent(0) &&
        (it_x - xvf.data()) == (it_y - yvf.data()))
    {
      size_t j = it_x - xvf.data();
      p0(i) = pf(j);
    }
  }

  return p0;
}

ViewVector_d Warmstart(const ViewVector_d& xv0, const ViewVector_d& yv0, const ViewVector_d& xvf,
    const ViewVector_d& yvf, const ViewVector_d& pf)
{
  return ViewVector_d(
      "void", 0);  // # TODO: to do this efficiently, we should store the indices that were removed,
                   // as they are removed, so in ContactSetPredictor or something idk
}
