#include "mirco_topologyutilities_kokkos.h"

#include <cmath>

ViewVector_d MIRCO::CreateMeshgrid(const int ngrid, const double GridSize)
{
  ViewVector_d meshgrid("meshgrid", ngrid);

  const double GridSize_2 = GridSize / 2;
  Kokkos::parallel_for(
      "Create meshgrid", ngrid,
      KOKKOS_LAMBDA(const int i) { meshgrid(i) = GridSize_2 + i * GridSize; });

  return meshgrid;
}

MIRCO::TopologyMaxAndMean MIRCO::ComputeMaxAndMean(ViewMatrix_d topology_d)
{
  const int n0 = topology_d.extent(0);
  const int n1 = topology_d.extent(1);

  // Note: these reductions could be merged into one, but that would require defining a custom
  // reducer struct type
  double zmax = -std::numeric_limits<double>::infinity();
  Kokkos::parallel_reduce(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {n0, n1}),
      KOKKOS_LAMBDA(int i, int j, double& update) {
        const double val = topology_d(i, j);
        if (val > update) update = val;
      },
      Kokkos::Max<double>(zmax));

  double zmean = 0.0;
  Kokkos::parallel_reduce(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {n0, n1}),
      KOKKOS_LAMBDA(int i, int j, double& update) { update += topology_d(i, j); }, zmean);
  zmean /= (n1 * n2);

  return TopologyMaxAndMean{zmax, zmean};
}
