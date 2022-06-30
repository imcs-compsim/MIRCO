#ifndef SRC_CONTACTPREDICTORS_H_
#define SRC_CONTACTPREDICTORS_H_

#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>

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
      std::vector<double> &b0, double zmax, double Delta, double w_el, std::vector<double> meshgrid,
      Epetra_SerialDenseMatrix topology);

  /**
   * @brief The aim of this function is to guess the set of nodes in contact among the nodes
   * predicted in the ContactSetPredictor function. It uses Warmstart to make an initial guess of
   * the nodes incontact in this iteration based on the previous iteration.
   *
   * @param flagwarm Warm-Starter flag
   * @param k Iteration number
   * @param n0 Number of nodes predicted to be in contact
   * @param nf Number of nodes in contact in the previous iteration
   * @param xv0 x-coordinates of the points in contact in the previous iteration.
   * @param yv0 y-coordinates of the points in contact in the previous iteration.
   * @param pf Contact force at (xvf,yvf) predicted in the previous iteration.
   * @param x0 contact forces at (xvf,yvf) predicted in the previous iteration but are a part of
   * currect predicted contact set.
   * @param b0 Indentation value of the half space at the predicted points of contact.
   * @param xvf x-coordinates of the points in contact in the previous iteration.
   * @param yvf y-coordinates of the points in contact in the previous iteration.
   */
  void InitialGuessPredictor(bool flagwarm, int k, int n0, int nf, std::vector<double> xv0,
      std::vector<double> yv0, std::vector<double> pf, Epetra_SerialDenseMatrix &x0,
      std::vector<double> &b0, std::vector<double> xvf, std::vector<double> yvf);
}  // namespace MIRCO

#endif  // SRC_CONTACTPREDICTORS_H_