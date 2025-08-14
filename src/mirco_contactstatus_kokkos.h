#ifndef SRC_CONTACTSTATUS_KOKKOS_H_
#define SRC_CONTACTSTATUS_KOKKOS_H_

#include "mirco_kokkostypes_kokkos.h"

namespace MIRCO
{
  /**
   * @brief The aim of this function is to compute to nodes in contact in the current iteration
   *
   * @param[out] xvf x-coordinates of the points in contact in the previous iteration.
   * @param[out] yvf y-coordinates of the points in contact in the previous iteration.
   * @param[out] pf Contact force at (xvf,yvf) predicted in the previous iteration.
   * @param[out] nf Number of nodes in contact in the previous iteration
   * @param[in] y Solution containing force
   * @param[in] xv0 x-coordinates of the points in contact in the previous iteration.
   * @param[in] yv0 y-coordinates of the points in contact in the previous iteration.
   */
  void ComputeContactNodes(ViewVector_d& xvf, ViewVector_d& yvf, ViewVector_d& pf_d,
      const int activeSetSize, const ViewVector_d p_d, const ViewVector_d xv0,
      const ViewVector_d yv0);

  /**
   * @brief Calculate the contact force and contact area for the
   * current iteration.
   *
   * @param[out] totalForce Total force
   * @param[out] contactArea Contact area
   * @param[in] pf Contact force vector.
   * @param[in] GridSize Grid size (length of each cell)
   * @param[in] LateralLength Lateral side of the surface [micrometers]
   * @param[in] PressureGreenFunFlag Flag to use Green function based on uniform pressure instead of
   * point force
   */
  void ComputeContactForceAndArea(double& totalForce, double& contactArea, const ViewVector_d pf_d,
      const double GridSize, const double LateralLength, const bool PressureGreenFunFlag);
}  // namespace MIRCO

#endif  // SRC_CONTACTSTATUS_KOKKOS_H_
