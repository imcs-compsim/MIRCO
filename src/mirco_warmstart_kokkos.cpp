#include "mirco_warmstart_kokkos.h"

#include <algorithm>

namespace MIRCO
{
  ViewVector_h Warmstart(const ViewVector_h& xv0, const ViewVector_h& yv0, const ViewVector_h& xvf,
      const ViewVector_h& yvf, const ViewVector_h& pf)
  {
    // TODO: If possible, convert this algorithm into running in a parallel kernel so that device
    // views can be used
    const auto n = static_cast<size_t>(xv0.extent(0));
    ViewVector_h p0("MIRCO::Warmstart(); p0", n);

    for (size_t i = 0; i < n; ++i)
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

}  // namespace MIRCO
