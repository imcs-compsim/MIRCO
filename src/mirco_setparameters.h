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
   * @param LateralLength Lateral side of the surface [micrometers]
   * @param nu1 Poisson's ratio of body 1
   * @param nu2 Poisson's ratio of body 2
   * @param CompositeYoungs The composite Young's modulus
   * @param alpha Correction factor for this problem. Depends on the resolution.
   * @param ElasticComplianceCorrection Elastic compliance correction
   * @param GridSize Grid size (length of each cell)
   * @param Tolerance Tolerance for the convergence of force.
   * @param Delta Far-field displacement (Gap).
   * @param TopologyFilePath Path of the input file containing the topology.
   * @param Resolution Resolution parameter
   * @param InitialTopologyStdDeviation Initial Standard deviation for the random-midpoint generator
   * [micrometers]
   * @param inputFileName The name of the input file containing all the parameters
   * @param RandomTopologyFlag Set `true` to use the Random-Midpoint Generator to generate the
   * topology. Set `false` to read topology from a file.
   * @param Hurst Hurst Exponent (Used in random mid-point generator)
   * @param RandomSeedFlag Set `true` to fix the seed to generate psuedo random topology to
   * reproduce results. Set `false` to use random seed.
   * @param RandomGeneratorSeed Set the value of seed for the random mid-point generator
   * @param WarmStartingFlag Set `true` for using the warm starter. It predicts the nodes coming
   * into contact in the next iteration and hence speeds up the computation.
   * @param MaxIteration Maximum number of iterations for the force to converge.
   * @param PressureGreenFunFlag Flag to use Green function based on uniform pressure instead of
   * point force
   * @param ExpectedPressure Expected pressure output for framework test
   * @param ExpectedPressureTolerance Tolerance to compare the output to the expected pressure
   */
  void SetParameters(double& E1, double& E2, double& LateralLength, double& nu1, double& nu2,
      double& CompositeYoungs, double& alpha, double& ElasticComplianceCorrection, double& GridSize,
      double& Tolerance, double& Delta, std::string& TopologyFilePath, int& Resolution,
      double& InitialTopologyStdDeviation, const std::string& inputFileName,
      bool& RandomTopologyFlag, double& Hurst, bool& RandomSeedFlag, int& RandomGeneratorSeed,
      bool& WarmStartingFlag, int& MaxIteration, bool& PressureGreenFunFlag,
      double& ExpectedPressure, double& ExpectedPressureTolerance);
}  // namespace MIRCO

#endif  // SRC_SETPARAMETERS_H_
