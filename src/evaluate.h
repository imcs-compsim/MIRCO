#ifndef SRC_EVALUATE_H_
#define SRC_EVALUATE_H_

#include <Epetra_SerialSymDenseMatrix.h>
#include <string>

namespace MIRCO
{
  /**
   * @brief Relate the far-field displacement with pressure
   *
   * @param pressure Pressure
   * @param Delta Far-field displacement (Gap)
   * @param lato Lateral side of the surface [micrometers]
   * @param delta Grid size
   * @param errf The error in the force vector
   * @param to1 Tolerance for the convergence of force
   * @param max_iter Maximum number of iterations for the force to converge
   * @param E Composite Young's modulus
   * @param flagwarm Warm-Starter Flag
   * @param k_el Elastic compliance correction
   * @param topology Topology matrix containing heights
   * @param zmax Maximum height
   * @param meshgrid Meshgrid vector
   */
  void Evaluate(double &pressure, double Delta, double lato, double delta, double errf, double to1,
      int max_iter, double E, bool flagwarm, double k_el, Epetra_SerialDenseMatrix topology,
      double zmax, std::vector<double> meshgrid);
}  // namespace MIRCO


#endif  // SRC_EVALUATE_H_
