#ifndef SRC_TOPOLOGYUTILITIES_H_
#define SRC_TOPOLOGYUTILITIES_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <cmath>
#include <vector>

#include <Kokkos_Core.hpp>
// tmp
#include "tmpHelpers/Timer.hpp"
#include "tmpHelpers/kokkosIntegration.hpp"

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
    double max_;  /*!< Maximum height of the topology */
    double mean_; /*!< Mean height of the topology */
  };

  /**
   * @brief Compute the maximum height and the mean height of the topology.
   *
   * @param topology Topology matrix containing heights
   */
  TopologyMaxAndMean ComputeMaxAndMean(const Teuchos::SerialDenseMatrix<int, double>& topology);
}  // namespace MIRCO

#endif  // SRC_TOPOLOGYUTILITIES_H_
