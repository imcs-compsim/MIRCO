#include "mirco_contactpredictors_kokkos.h"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <vector>

#include "mirco_warmstart_kokkos.h"

void MIRCO::ContactSetPredictor(int &n0, ViewVector_d &xv0, ViewVector_d &yv0, ViewVector_d &b0,
    double zmax, double Delta, double w_el, const ViewVector_d meshgrid,
    const ViewMatrix_d topology)
{
  // # TODO:m see if this parallization is even needed or just makes it less efficient

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

  Kokkos::View<int*, Device_Default_t> row("row", n0);
  Kokkos::View<int*, Device_Default_t> col("col", n0);

  //# might need to do int* but just of size 1, because we need a ptr on device
  ///Kokkos::View<int, Device_Default_t> counter("counter");
  Kokkos::View<int*, Device_Default_t> counter("counter", 1);
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
        
        // b0 = \overbar{u} + w_el
        b0(i) = Delta + w_el - (zmax - topology(ri, ci));
      });
}
