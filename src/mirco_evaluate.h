#ifndef SRC_EVALUATE_H_
#define SRC_EVALUATE_H_

#include <Teuchos_SerialDenseMatrix.hpp>
#include <string>

namespace MIRCO
{
  /**
   * @brief Relate the far-field displacement with pressure
   *
   * @param Delta Far-field displacement (Gap)
   * @param LateralLength Lateral side of the surface [micrometers]
   * @param GridSize Grid size
   * @param Tolerance Tolerance for the convergence of force
   * @param MaxIteration Maximum number of iterations for the force to converge
   * @param CompositeYoungs Composite Young's modulus
   * @param WarmStartingFlag Warm-Starter Flag
   * @param ElasticComplianceCorrection Elastic compliance correction
   * @param topology Topology matrix containing heights
   * @param zmax Maximum height
   * @param meshgrid Meshgrid vector
   * @param xvf x-Coordinate vector of points that are in contact, calculated in each iteration
   * @param yvf y-Coordinate vector of points that are in contact, calculated in each iteration
   * @param pf Contact force vector at (xvf,yvf) calculated in each iteration.
   * @param nf Number of contact points
   */
  void Evaluate(double Delta, double LateralLength, double GridSize, double Tolerance,
      int MaxIteration, double CompositeYoungs, bool WarmStartingFlag,
      double ElasticComplianceCorrection, Teuchos::SerialDenseMatrix<int, double>& topology,
      double zmax, std::vector<double>& meshgrid, std::vector<double>& xvf,
      std::vector<double>& yvf, std::vector<double>& pf, int& n);
}  // namespace MIRCO


#endif  // SRC_EVALUATE_H_
