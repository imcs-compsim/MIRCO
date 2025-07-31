#ifndef SRC_CONTACTSTATUS_H_
#define SRC_CONTACTSTATUS_H_

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
  void ComputeContactNodes(ViewVector_d& xvf, ViewVector_d& yvf,
    ViewVector_d& pf_d, const int activeSetSize, const ViewVector_d p_d,
    const ViewVector_d xv0, const ViewVector_d yv0);

  /**
   * @brief The aim of this function is to calulate the contact force and contact area for the
   * current iteration. Note that the force vector
   *
   * @param[out] totalForce Total force
   * @param[out] contactArea Contact area
   * @param nf Number of nodes in contact in the previous iteration
   * @param pf Contact force vector predicted in the previous iteration.
   * @param GridSize Grid size (length of each cell)
   * @param LateralLength Lateral side of the surface [micrometers]
   * @param PressureGreenFunFlag Flag to use Green function based on uniform pressure instead of
   * point force
   */
  void ComputeContactForceAndArea(double& totalForce, double& contactArea, const ViewVector_d pf_d, const double GridSize, const double LateralLength, const bool PressureGreenFunFlag);
}  // namespace MIRCO

#endif  // SRC_CONTACTSTATUS_H_
