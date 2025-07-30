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
   * @param[out] yv0 y-coordinates of the points in contact in the previous iteration. //# TODO:g
   * consider switching the order of in and out params and also keep it consistent everywhere
   * @param[out] b0 Indentation value of the half space at the predicted points of contact.
   * @param[in] zmax Maximum height of the topology
   * @param[in] Delta Far-field displacement (Gap)
   * @param[in] w_el Elastic correction
   * @param[in] meshgrid Meshgrid
   * @param[in] topology Topology matrix containing heights
   */
  void ContactSetPredictor(int &n0, ViewVector_d &xv0, ViewVector_d &yv0, ViewVector_d &b0,
      double zmax, double Delta, double w_el, const ViewVector_d &meshgrid,
      const ViewMatrix_d &topology);

  /**
   * @brief The aim of this function is to guess the set of nodes in contact among the nodes
   * predicted in the ContactSetPredictor function. It uses Warmstart to make an initial guess of
   * the nodes incontact in this iteration based on the previous iteration.
   *
   * @param[in] WarmStartingFlag Warm-Starter flag
   * @param[in] k Iteration number
   * @param[in] n0 Number of nodes predicted to be in contact
   * @param[in] xv0 x-coordinates of the points in contact in the previous[in] iteration.
   * @param[in] yv0 y-coordinates of the points in contact in the previous iteration.
   * @param[in] pf Contact force at (xvf,yvf) predicted in the previous iteration.
   * @param[in] x0 contact forces at (xvf,yvf) predicted in the previous iteration but are a part of
   * currect predicted contact set.
   * @param[in] b0 Indentation value of the half space at the predicted points of contact.
   * @param[in] xvf Coordinates of the points predicted to be in contact. //#
   * @param[in] yvf Coordinates of the points predicted to be in contact. //#
   *
   * @return p0 vector of contact forces predicted in the previous iteration but are a part of
   * currect predicted contact set.
   */
  ViewVector_h InitialGuessPredictor(const bool WarmStartingFlag, const int k, const int n0, const ViewVector_h &xv0,
      const ViewVector_h &yv0, const ViewVector_h &pf, const ViewVector_h &b0,
      const ViewVector_h &xvf, const ViewVector_h &yvf);
}  // namespace MIRCO

#endif  // SRC_CONTACTPREDICTORS_KOKKOS_H_
