#ifndef SRC_EVALUATE_H_
#define SRC_EVALUATE_H_

#include <Kokkos_Core.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>

#include "mirco_inputparameters.h"
// tmp
#include "tmpHelpers/Timer.hpp"
#include "tmpHelpers/kokkosIntegration.hpp"

namespace MIRCO
{
  /**
   * @brief Relate the far-field displacement with pressure
   *
   * @param pressure Pressure
   * @param Delta Far-field displacement (Gap)
   * @param LateralLength Lateral side of the surface [micrometers]
   * @param GridSize Grid size (length of each cell)
   * @param Tolerance Tolerance for the convergence of force
   * @param MaxIteration Maximum number of iterations for the force to converge
   * @param CompositeYoungs Composite Young's modulus
   * @param WarmStartingFlag Warm-Starter Flag
   * @param ElasticComplianceCorrection Elastic compliance correction
   * @param topology Topology matrix containing heights
   * @param zmax Maximum height
   * @param meshgrid Meshgrid vector
   * @param PressureGreenFunFlag Flag to use Green function based on uniform pressure instead of
   * point force
   */
  void Evaluate(double& pressure, const double Delta, const double LateralLength,
      const double GridSize, const double Tolerance, const int MaxIteration,
      const double CompositeYoungs, const bool WarmStartingFlag,
      const double ElasticComplianceCorrection,
      const Teuchos::SerialDenseMatrix<int, double>& topology, const double zmax,
      const std::vector<double>& meshgrid, const bool PressureGreenFunFlag);

  /**
   * @brief Relate the far-field displacement with pressure, taking the parameters from a
   * MIRCO::InputParameters object
   *
   * @param inputParams Object which holds the input parameters
   * @param pressure Pressure
   * @param zmax Maximum height
   * @param meshgrid Meshgrid vector
   */
  inline void Evaluate(const MIRCO::InputParameters& inputParams, double& pressure,
      const double zmax, const std::vector<double>& meshgrid)
  {
    Evaluate(pressure, inputParams.delta_, inputParams.lateral_length_, inputParams.grid_size_,
        inputParams.tolerance_, inputParams.max_iteration_, inputParams.composite_youngs_,
        inputParams.warm_starting_flag_, inputParams.elastic_compliance_correction_,
        *(inputParams.topology_), zmax, meshgrid, inputParams.pressure_green_funct_flag_);
  }
}  // namespace MIRCO


#endif  // SRC_EVALUATE_H_
