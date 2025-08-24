#include "mirco_contactstatus.h"

#include <cmath>

namespace MIRCO
{
  void ComputeContactForceAndArea(double& totalForce, double& contactArea, const ViewVector_d pf,
      const double GridSize, const double LateralLength, const bool PressureGreenFunFlag)
  {
    totalForce = 0;

    const int activeSetSize = pf.extent(0);
    const double GridSize2 = GridSize * GridSize;

    if (PressureGreenFunFlag)
      Kokkos::parallel_reduce(
          activeSetSize,
          KOKKOS_LAMBDA(const int i, double& local_sum) { local_sum += pf(i) * GridSize2; },
          totalForce);
    else
      Kokkos::parallel_reduce(
          activeSetSize, KOKKOS_LAMBDA(const int i, double& local_sum) { local_sum += pf(i); },
          totalForce);

    contactArea = activeSetSize * (GridSize2 / LateralLength * LateralLength);
  }

}  // namespace MIRCO
