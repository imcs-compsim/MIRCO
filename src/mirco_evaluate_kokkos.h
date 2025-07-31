#ifndef SRC_EVALUATE_KOKKOS_H_
#define SRC_EVALUATE_KOKKOS_H_

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
   * @param[in] meshgrid Meshgrid vector
   * @param[in] PressureGreenFunFlag Flag to use Green function based on uniform pressure instead of
   * point force
   *
   * @return Mean pressure
   */
  double Evaluate(const double Delta, const double LateralLength,
      const double GridSize, const double Tolerance, const int MaxIteration,
      const double CompositeYoungs, const bool WarmStartingFlag,
      const double ElasticComplianceCorrection,
      const Teuchos::SerialDenseMatrix<int, double>& topology, const double zmax,
      const ViewVector_d meshgrid_d, const bool PressureGreenFunFlag);

  /**
   * @brief Relate the far-field displacement with pressure, taking the parameters from a
   * MIRCO::InputParameters object
   *
   * @param[in] inputParams Object which holds the input parameters
   * @param[in] zmax Maximum height
   * @param[in] meshgrid Meshgrid vector
   *
   * @return Mean pressure
   */
  inline double Evaluate(const MIRCO::InputParameters& inputParams,
      const double zmax, const ViewVector_d meshgrid_d)
  {
    return Evaluate(inputParams.delta_, inputParams.lateral_length_, inputParams.grid_size_,
        inputParams.tolerance_, inputParams.max_iteration_, inputParams.composite_youngs_,
        inputParams.warm_starting_flag_, inputParams.elastic_compliance_correction_,
        *(inputParams.topology_), zmax, meshgrid_d, inputParams.pressure_green_funct_flag_);
  }
}  // namespace MIRCO

#endif  // SRC_EVALUATE_KOKKOS_H_
