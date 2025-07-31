#ifndef SRC_CONTACTPREDICTORS_KOKKOS_H_
#define SRC_CONTACTPREDICTORS_KOKKOS_H_

#include <Kokkos_Core.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <vector>
// tmp
#include "tmpHelpers/Timer.hpp"
#include "tmpHelpers/kokkosIntegration.hpp"

namespace MIRCO
{
  /**
   * @brief The aim of this function is to determine all the points, for which the gap is bigger
   * than the displacement of the rigid indenter, cannot be in contact and thus are not checked in
   * nonlinear solve
   *
   * @param[out] n0 Number of nodes predicted to be in contact
   * @param[out] xv0 x-coordinates of the points in contact in the previous iteration.
   * @param[out] yv0 y-coordinates of the points in contact in the previous iteration.
   * @param[out] b0 Indentation value of the half space at the predicted points of contact.
   * @param[in] zmax Maximum height of the topology
   * @param[in] Delta Far-field displacement (Gap)
   * @param[in] w_el Elastic correction
   * @param[in] meshgrid Meshgrid
   * @param[in] topology Topology matrix containing heights
   */
  void ContactSetPredictor(int &n0, ViewVector_d &xv0, ViewVector_d &yv0, ViewVector_d &b0,
      double zmax, double Delta, double w_el, const ViewVector_d meshgrid,
      const ViewMatrix_d topology);
}  // namespace MIRCO

#endif  // SRC_CONTACTPREDICTORS_KOKKOS_H_
