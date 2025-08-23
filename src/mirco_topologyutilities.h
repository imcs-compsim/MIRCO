#ifndef SRC_TOPOLOGYUTILITIES_KOKKOS_H_
#define SRC_TOPOLOGYUTILITIES_KOKKOS_H_

#include "mirco_kokkostypes_kokkos.h"

namespace MIRCO
{
  /**
   * @brief Create a Meshgrid vector of the surface
   *
   * Creates a meshgrid with coordinates. Since the vector is identical in both directions,
   * therefore only one is created.
   *
   * @param meshgrid Meshgrid vector
   * @param ngrid Number of grid points in one direction
   * @param GridSize Grid size (length of each cell)
   */
  ViewVector_d CreateMeshgrid(const int ngrid, const double GridSize);

  /**
   * @brief Store maximum and mean height of the topology
   *
   */
  struct TopologyMaxAndMean
  {
    double max;  /*!< Maximum height of the topology */
    double mean; /*!< Mean height of the topology */
  };

  /**
   * @brief Compute the maximum height and the mean height of the topology.
   *
   * @param topology Topology matrix containing heights
   */
  TopologyMaxAndMean ComputeMaxAndMean(ViewMatrix_d topology_d);
}  // namespace MIRCO

#endif  // SRC_TOPOLOGYUTILITIES_KOKKOS_H_
