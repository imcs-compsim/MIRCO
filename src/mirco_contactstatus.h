#ifndef SRC_CONTACTSTATUS_H_
#define SRC_CONTACTSTATUS_H_

#include <Teuchos_SerialDenseMatrix.hpp>
#include <cmath>
#include <vector>

namespace MIRCO
{
  /**
   * @brief The aim of this function is to compute to nodes in contact in the current iteration
   *
   * @param xvf x-coordinates of the points in contact in the previous iteration.
   * @param yvf y-coordinates of the points in contact in the previous iteration.
   * @param pf Contact force at (xvf,yvf) predicted in the previous iteration.
   * @param nf Number of nodes in contact in the previous iteration
   * @param y Solution containing force
   * @param xv0 x-coordinates of the points in contact in the previous iteration.
   * @param yv0 y-coordinates of the points in contact in the previous iteration.
   */
  void ComputeContactNodes(std::vector<double> &xvf, std::vector<double> &yvf,
      std::vector<double> &pf, int &nf, Teuchos::SerialDenseMatrix<int, double> y,
      std::vector<double> xv0, std::vector<double> yv0);

  /**
   * @brief The aim of this function is to calulate the contact force and contact area for the
   * current iteration
   *
   * @param force0 Force vector; Each element contating contact force calculated at every iteraion
   * @param area0 Force vector; Each element contating contact area calculated at every iteraion
   * @param w_el Elastic correction
   * @param nf Number of nodes in contact in the previous iteration
   * @param pf Contact force at (xvf,yvf) predicted in the previous iteration.
   * @param k Iteration number
   * @param GridSize Grid size
   * @param LateralLength Lateral side of the surface [micrometers]
   * @param ElasticComplianceCorrection Elastic compliance correction
   */
  void ComputeContactForceAndArea(std::vector<double> &force0, std::vector<double> &area0,
      double &w_el, int nf, std::vector<double> pf, int k, double GridSize, double LateralLength,
      double ElasticComplianceCorrection);
}  // namespace MIRCO

#endif  // SRC_CONTACTSTATUS_H_