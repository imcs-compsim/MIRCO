#include "mirco_contactstatus_kokkos.h"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <cmath>
#include <vector>

void MIRCO::ComputeContactNodes(ViewVector_d& xvf, ViewVector_d& yvf,
    ViewVector_d& pf_d, const int activeSetSize, const ViewVector_d p_d,
    const ViewVector_d xv0, const ViewVector_d yv0)
{
  xvf = ViewVector_d("new xvf", activeSetSize);
  yvf = ViewVector_d("new xvf", activeSetSize);
  pf_d = ViewVector_d("new xvf", activeSetSize);
  int n0 = p_d.extent(0);
  // @} Parallelizing this slows down program, so removed it.
  
  Kokkos::View<int*, Device_Default_t> counter("counter", 1);
  Kokkos::deep_copy(counter, 0);

  Kokkos::parallel_for("", Kokkos::RangePolicy<ExecSpace_Default_t>(0, n0),
    KOKKOS_LAMBDA(const int i) {
      if(p_d(i) != 0) {
        int pos = Kokkos::atomic_fetch_add(&counter(0), 1);
        xvf(pos) = xv0(i);
        yvf(pos) = yv0(i);
        pf_d(pos) = p_d(i);
      }
    });
}

void MIRCO::ComputeContactForceAndArea(double& totalForce, double& contactArea, const ViewVector_d pf_d, const double GridSize, const double LateralLength, const bool PressureGreenFunFlag)
{
  const int activeSetSize = pf_d.extent(0);
  totalForce = 0;
  
  const double GridSize2 = GridSize*GridSize;
  
  if(PressureGreenFunFlag)
    Kokkos::parallel_reduce("", activeSetSize,
      KOKKOS_LAMBDA(int i, double& local_sum) {
        local_sum += pf_d(i) * GridSize2;
      },
      totalForce);
  else
    Kokkos::parallel_reduce("", activeSetSize,
      KOKKOS_LAMBDA(int i, double& local_sum) {
        local_sum += pf_d(i);
      },
      totalForce);
      
      
  contactArea = activeSetSize * (GridSize2 / LateralLength*LateralLength);
}
