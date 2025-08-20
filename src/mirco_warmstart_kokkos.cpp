#include "mirco_warmstart_kokkos.h"

#include <algorithm>

namespace MIRCO
{
  ViewVector_d Warmstart(
      const ViewVector_d& activeSet0_d, const ViewVector_d& activeSetf_d, const ViewVector_d& pf_d)
  {
    const int n0 = activeSet0_d.extent(0);
    const int nf = activeSetf_d.extent(0);
    ViewVector_d p0_d("Warmstart(); p0", n0);

    Kokkos::parallel_for(
        n0, KOKKOS_LAMBDA(const int i) {
          const int a0_i = activeSet0_d(i);
          for (int j = 0; j < nf; ++j)
          {
            if (activeSetf_d(j) == a0_i)
            {
              p0_d(i) = pf_d(j);
              break;
            }
          }
        });

    return p0_d;
  }

}  // namespace MIRCO
