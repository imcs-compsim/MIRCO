#ifndef SRC_CONTACTSTATUS_KOKKOS_H_
#define SRC_CONTACTSTATUS_KOKKOS_H_

#include "mirco_kokkostypes_kokkos.h"

namespace MIRCO
{
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
