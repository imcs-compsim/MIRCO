#include "mirco_contactpredictors_kokkos.h"

namespace MIRCO
{
  void ContactSetPredictor(int &n0, ViewVector_d &xv0, ViewVector_d &yv0, ViewVector_d &b0,
      double zmax, double Delta, double w_el, const ViewVector_d meshgrid,
      const ViewMatrix_d topology)
  {
    n0 = 0;

    double value = zmax - Delta - w_el;

    int N = topology.extent(0);

    Kokkos::parallel_reduce(
        "ContactSetPredictor0", N * N,
        KOKKOS_LAMBDA(int idx, int &local_sum) {
          int i = idx / N;
          int j = idx % N;
          if (topology(i, j) >= value) local_sum++;
        },
        n0);

    ViewVectorInt_d row("row", n0);
    ViewVectorInt_d col("col", n0);

    ViewVectorInt_d counter("counter", 1);
    Kokkos::deep_copy(counter, 0);

    Kokkos::parallel_for(
        "ContactSetPredictor1", N * N, KOKKOS_LAMBDA(int idx) {
          int i = idx / N;
          int j = idx % N;
          if (topology(i, j) >= value)
          {
            int pos = Kokkos::atomic_fetch_add(&counter(0), 1);
            row(pos) = i;
            col(pos) = j;
          }
        });

    xv0 = ViewVector_d("xv0", n0);
    yv0 = ViewVector_d("yv0", n0);
    b0 = ViewVector_d("b0", n0);

    Kokkos::parallel_for(
        "ContactSetPredictor2", n0, KOKKOS_LAMBDA(int i) {
          int ci = col(i);
          int ri = row(i);

          xv0(i) = meshgrid(ci);
          yv0(i) = meshgrid(ri);

          // Note: b0 = \overbar{u} + w_el
          b0(i) = Delta + w_el - (zmax - topology(ri, ci));
        });
  }

}  // namespace MIRCO
