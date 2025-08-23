#include "mirco_contactpredictors_kokkos.h"

namespace MIRCO
{
  void ContactSetPredictor(ViewVectorInt_d& activeSet0_d, ViewVector_d& xv0_d, ViewVector_d& yv0_d,
      ViewVector_d& b0, double zmax, double Delta, double w_el, const ViewMatrix_d topology_d,
      const ViewVector_d meshgrid_d)
  {
    const int N = topology_d.extent(0);

    const double critValue = zmax - Delta - w_el;
    int n0 = 0;
    Kokkos::parallel_reduce(
        N * N,
        KOKKOS_LAMBDA(const int a, int& local_sum) {
          if (topology_d(a / N, a % N) >= critValue) local_sum++;
        },
        n0);

    activeSet0_d = ViewVectorInt_d("activeSet0_d", n0);
    xv0_d = ViewVector_d("xv0_d", n0);
    yv0_d = ViewVector_d("yv0_d", n0);
    b0 = ViewVector_d("b0", n0);

    ViewScalarInt_d counter("ContactSetPredictor(); counter");
    Kokkos::deep_copy(counter, 0);
    const double Dwz = Delta + w_el - zmax;
    Kokkos::parallel_for(
        N * N, KOKKOS_LAMBDA(const int a) {
          const int i = a / N;
          const int j = a % N;
          const double topology_a = topology_d(i, j);
          if (topology_a >= critValue)
          {
            const int aa = Kokkos::atomic_fetch_add(&counter(), 1);
            activeSet0_d(aa) = a;
            xv0_d(aa) = meshgrid_d(i);
            yv0_d(aa) = meshgrid_d(j);
            // Note: b0 = \overbar{u} + w_el = Delta + w_el - (zmax - topology_a);
            b0(aa) = Dwz + topology_a;
          }
        });
  }

}  // namespace MIRCO
