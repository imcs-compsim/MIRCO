#ifndef SRC_INPUTPARAMETERS_KOKKOS_H_
#define SRC_INPUTPARAMETERS_KOKKOS_H_

#include <string>

#include "mirco_kokkostypes_kokkos.h"

namespace MIRCO
{
  /**
   * @brief This struct stores the input parameters and topology
   *
   */
  struct InputParameters
  {
    /**
     * @brief Constructor which sets the necessary member variable parameters from an input (.xml)
     * file and creates the topology
     *
     * @param inputFileName Input file w.r.t. the calling directory
     */
    InputParameters(const std::string& inputFileName);

    /**
     * @brief Constructor which sets the necessary member variable parameters without an input
     * (.xml) file and creates the topology using the random midpoint generator
     *
     * @param E1 Young's modulus of body 1
     * @param E2 Young's modulus of body 2
     * @param nu1 Poisson's ratio of body 1
     * @param nu2 Poisson's ratio of body 2
     * @param Tolerance Tolerance for the convergence of force.
     * @param Delta Far-field displacement (Gap).
     * @param LateralLength Lateral side of the surface [micrometers]
     * @param Resolution Resolution parameter
     * @param InitialTopologyStdDeviation Initial Standard deviation for the random-midpoint
     * generator [micrometers]
     * @param Hurst Hurst Exponent (Used in random mid-point generator)
     * @param RandomSeedFlag Set `true` to fix the seed to generate psuedo random topology to
     * reproduce results. Set `false` to use random seed.
     * @param RandomGeneratorSeed Set the value of seed for the random mid-point generator
     * @param MaxIteration Maximum number of iterations for the force to converge.
     * @param WarmStartingFlag Set `true` for using the warm starter. It predicts the nodes coming
     * into contact in the next iteration and hence speeds up the computation.
     * @param PressureGreenFunFlag Flag to use Green function based on uniform pressure instead of
     * point force
     */
    InputParameters(double E1, double E2, double nu1, double nu2, double Tolerance, double Delta,
        double LateralLength, int Resolution, double InitialTopologyStdDeviation, double Hurst,
        bool RandomSeedFlag, int RandomGeneratorSeed, int MaxIteration, bool WarmStartingFlag,
        bool PressureGreenFunFlag);

    /**
     * @brief Constructor which sets the necessary member variable parameters without an input
     * (.xml) file and creates the topology from a specified topology (.dat) file
     *
     * @param E1 Young's modulus of body 1
     * @param E2 Young's modulus of body 2
     * @param nu1 Poisson's ratio of body 1
     * @param nu2 Poisson's ratio of body 2
     * @param Tolerance Tolerance for the convergence of force.
     * @param Delta Far-field displacement (Gap).
     * @param LateralLength Lateral side of the surface [micrometers]
     * @param TopologyFilePath Path of the input file containing the topology.
     * @param MaxIteration Maximum number of iterations for the force to converge.
     * @param WarmStartingFlag Set `true` for using the warm starter. It predicts the nodes coming
     * into contact in the next iteration and hence speeds up the computation.
     * @param PressureGreenFunFlag Flag to use Green function based on uniform pressure instead of
     * point force
     */
    InputParameters(double E1, double E2, double nu1, double nu2, double Tolerance, double Delta,
        double LateralLength, const std::string& TopologyFilePath, int MaxIteration,
        bool WarmStartingFlag, bool PressureGreenFunFlag);

    int N = 0;
    double composite_youngs = 0.0, elastic_compliance_correction = 0.0, shape_factor = 0.0,
           tolerance = 0.0, delta = 0.0, lateral_length = 0.0, grid_size = 0.0;
    int max_iteration = 0;
    bool warm_starting_flag = false;
    bool pressure_green_funct_flag = false;
    // Note: topology_d is a lightweight handle, similar to std::shared_ptr. This struct does not
    // own topology_d.
    ViewMatrix_d topology_d;
  };
}  // namespace MIRCO

#endif  // SRC_INPUTPARAMETERS_KOKKOS_H_
