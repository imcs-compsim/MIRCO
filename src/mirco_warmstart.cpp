#include "mirco_warmstart.h"

namespace MIRCO
{
  ViewVector_d Warmstart(
      const ViewVectorInt_d& activeSet0, const ViewVectorInt_d& activeSetf, const ViewVector_d& pf)
  {
    const int n0 = activeSet0.extent(0);
    const int nf = activeSetf.extent(0);
    ViewVector_d p0("p0", n0);

    Kokkos::parallel_for(
        n0, KOKKOS_LAMBDA(const int i) {
          const int a0_i = activeSet0(i);
          for (int j = 0; j < nf; ++j)
          {
            if (activeSetf(j) == a0_i)
            {
              p0(i) = pf(j);
              break;
            }
          }
        });

    return p0;
  }

}  // namespace MIRCO
