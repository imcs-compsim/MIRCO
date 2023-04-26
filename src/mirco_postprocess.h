#ifndef SRC_POSTPROCESS_H_
#define SRC_POSTPROCESS_H_

#include <Teuchos_SerialDenseMatrix.hpp>
#include <string>

namespace MIRCO
{
  /**
   * @brief Postprocessing of contact force and contact area
   *
   * @param xvf x-Coordinate vector of points that are in contact, calculated in each iteration
   * @param yvf y-Coordinate vector of points that are in contact, calculated in each iteration
   * @param pf Contact force vector at (xvf,yvf) calculated in each iteration.
   * @param nf Number of contact points
   * @param GridSize Grid size
   * @param ngrid Number of grids
   * @param meshgrid Meshgrid vector
   * @param LateralLength Lateral side of the surface [micrometers]
   * @param elapsedTimeString Number of contact points
   * @param outputFileName Output file to store results
   */
  void PostProcess(std::vector<double> xvf, std::vector<double> yvf, std::vector<double> pf,
      int nf, double GridSize, int ngrid, std::vector<double>& meshgrid, double LateralLength,
      double elapsedTime, std::string outputFileName);
}  // namespace MIRCO


#endif  // SRC_POSTPROCESS_H_
