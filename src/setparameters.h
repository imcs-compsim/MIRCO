#ifndef SRC_SETPARAMETERS_H_
#define SRC_SETPARAMETERS_H_

#include <string>
namespace MIRCO
{
  /**
   * @brief Set the Parameters
   *
   * The aim of this fuction is to read the input file containing the simulation
   * specific, material and geometrical parameters.
   *
   * @param E1 Young's modulus of body 1
   * @param E2 Young's modulus of body 2
   * @param lato Lateral side of the surface [micrometers]
   * @param nu1 Poisson's ratio of body 1
   * @param nu2 Poisson's ratio of body 2
   * @param G1 Shear Modulus of body 1
   * @param G2 Shear Modulus of body 1
   * @param E The composite Young's modulus
   * @param alpha Correction factor for this problem. Depends on the resolution.
   * @param k_el Elastic compliance correction
   * @param delta Grid size
   * @param nnodi Total number of nodes
   * @param errf The error in the force vector. Initialise it with a high value.
   * @param tol Tolerance for the convergence of force.
   * @param Delta Far-field displacement (Gap).
   * @param zfilePath Path of the input file containing the topology.
   * @param resolution Resolution parameter
   * @param user_zmax Maximum height of the topology
   * @param inputFileName The name of the input file containing all the parameters
   * @param rmg_flag Set `true` to use the Random-Midpoint Generator to generate the topology. Set
   * `false` to read topology from a file.
   * @param Hurst Hurst Exponent (Used in random mid-point generator)
   * @param rand_seed_flag Set `true` to fix the seed to generate psuedo random topology to
   * reproduce results. Set `false` to use random seed.
   * @param rmg_seed Set the value of seed for the random mid-point generator
   * @param flagwarm Set `true` for using the warm starter. It predicts the nodes coming into
   * contact in the next iteration and hence speeds up the computation.
   * @param max_iter Maximum number of iterations for the force to converge.
   */
  void SetParameters(double& E1, double& E2, double& lato, double& nu1, double& nu2, double& G1,
      double& G2, double& E, double& alpha, double& k_el, double& delta, double& nnodi,
      double& errf, double& tol, double& Delta, std::string& zfilePath, int& resolution,
      double& user_zmax, const std::string& inputFileName, bool& rmg_flag, double& Hurst,
      bool& rand_seed_flag, int& rmg_seed, bool& flagwarm, int& max_iter);
}  // namespace MIRCO

#endif  // SRC_SETPARAMETERS_H_
