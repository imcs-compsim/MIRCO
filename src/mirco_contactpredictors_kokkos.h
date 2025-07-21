#ifndef SRC_CONTACTPREDICTORS_H_
#define SRC_CONTACTPREDICTORS_H_

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
   * @param n0 Number of nodes predicted to be in contact
   * @param xv0 x-coordinates of the points in contact in the previous iteration.
   * @param yv0 y-coordinates of the points in contact in the previous iteration.
   * @param b0 Indentation value of the half space at the predicted points of contact.
   * @param zmax Maximum height of the topology
   * @param Delta Far-field displacement (Gap)
   * @param w_el Elastic correction
   * @param meshgrid Meshgrid
   * @param topology Topology matrix containing heights
   */
  void ContactSetPredictor(int &n0, std::vector<double> &xv0, std::vector<double> &yv0,
      std::vector<double> &b0, double zmax, double Delta, double w_el,
      const std::vector<double> &meshgrid, const Teuchos::SerialDenseMatrix<int, double> &topology);

  /**
   * @brief The aim of this function is to guess the set of nodes in contact among the nodes
   * predicted in the ContactSetPredictor function. It uses Warmstart to make an initial guess of
   * the nodes incontact in this iteration based on the previous iteration.
   *
   * @param WarmStartingFlag Warm-Starter flag
   * @param k Iteration number
   * @param n0 Number of nodes predicted to be in contact
   * @param xv0 x-coordinates of the points in contact in the previous iteration.
   * @param yv0 y-coordinates of the points in contact in the previous iteration.
   * @param pf Contact force at (xvf,yvf) predicted in the previous iteration.
   * @param x0 contact forces at (xvf,yvf) predicted in the previous iteration but are a part of
   * currect predicted contact set.
   * @param b0 Indentation value of the half space at the predicted points of contact.
   * @param xvf Coordinates of the points predicted to be in contact. //#
   * @param yvf Coordinates of the points predicted to be in contact. //#
   */
  void InitialGuessPredictor(bool WarmStartingFlag, int k, int n0, const std::vector<double> &xv0,
      const std::vector<double> &yv0, const std::vector<double> &pf,
      Teuchos::SerialDenseMatrix<int, double> &x0, const std::vector<double> &b0,
      const std::vector<double> &xvf, const std::vector<double> &yvf);
}  // namespace MIRCO

#endif  // SRC_CONTACTPREDICTORS_H_
