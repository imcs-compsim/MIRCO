#ifndef SRC_CONTACTPREDICTORS_H_
#define SRC_CONTACTPREDICTORS_H_

#include "mirco_kokkostypes.h"

namespace MIRCO
{
  /**
   * @brief The aim of this function is to determine all the points, for which the gap is bigger
   * than the displacement of the rigid indenter, cannot be in contact and thus are not checked in
   * nonlinear solve
   *
   * @param[out] activeSet0 Points predicted to be in contact in the current iteration
   * @param[out] xv0 x-coordinate of the points predicted to be in contact in the current
   * iteration
   * @param[out] yv0 y-coordinate of the points predicted to be in contact in the current
   * iteration
   * @param[out] b0 Indentation value of the half space at the predicted points of contact
   * @param[in] zmax Maximum height of the topology
   * @param[in] Delta Far-field displacement (Gap)
   * @param[in] w_el Elastic correction
   * @param[in] topology Topology matrix containing heights
   * @param[in] meshgrid Meshgrid (coordinates in one direction)
   */
  void ContactSetPredictor(ViewVectorInt_d& activeSet0, ViewVector_d& xv0, ViewVector_d& yv0,
      ViewVector_d& b0, double zmax, double Delta, double w_el, const ViewMatrix_d topology,
      const ViewVector_d meshgrid);
}  // namespace MIRCO

#endif  // SRC_CONTACTPREDICTORS_H_
