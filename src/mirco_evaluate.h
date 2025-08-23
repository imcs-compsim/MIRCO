#ifndef SRC_EVALUATE_KOKKOS_H_
#define SRC_EVALUATE_KOKKOS_H_

#include "mirco_inputparameters_kokkos.h"
#include "mirco_kokkostypes_kokkos.h"

namespace MIRCO
{
  /**
   * @brief Relate the far-field displacement with pressure
   *
   * @param[out] pressure Mean pressure
   * @param[in] Delta Far-field displacement (Gap)
   * @param[in] LateralLength Lateral side of the surface [micrometers]
   * @param[in] GridSize Grid size (length of each cell)
   * @param[in] Tolerance Tolerance for the convergence of force
   * @param[in] MaxIteration Maximum number of iterations for the force to converge
   * @param[in] CompositeYoungs Composite Young's modulus
   * @param[in] WarmStartingFlag Warm-Starter Flag
   * @param[in] ElasticComplianceCorrection Elastic compliance correction
   * @param[in] topology Topology matrix containing heights
   * @param[in] zmax Maximum height
   * @param[in] meshgrid_d Meshgrid vector
   * @param[in] PressureGreenFunFlag Flag to use Green function based on uniform pressure instead of
   * point force
   */
  void Evaluate(double& pressure, const double Delta, const double LateralLength,
      const double GridSize, const double Tolerance, const int MaxIteration,
      const double CompositeYoungs, const bool WarmStartingFlag,
      const double ElasticComplianceCorrection, const ViewMatrix_d topology_d, const double zmax,
      const ViewVector_d meshgrid_d, const bool PressureGreenFunFlag);

  /**
   * @brief Relate the far-field displacement with pressure, taking the parameters from an
   * InputParameters object
   *
   * @param[out] pressure Mean pressure
   * @param[in] inputParams Object which holds the input parameters
   * @param[in] zmax Maximum height
   * @param[in] meshgrid_d Meshgrid vector
   */
  inline void Evaluate(double& pressure, const InputParameters& inputParams, const double zmax,
      const ViewVector_d meshgrid_d)
  {
    Evaluate(pressure, inputParams.delta, inputParams.lateral_length, inputParams.grid_size,
        inputParams.tolerance, inputParams.max_iteration, inputParams.composite_youngs,
        inputParams.warm_starting_flag, inputParams.elastic_compliance_correction,
        inputParams.topology_d, zmax, meshgrid_d, inputParams.pressure_green_funct_flag);
  }
}  // namespace MIRCO

#endif  // SRC_EVALUATE_KOKKOS_H_
