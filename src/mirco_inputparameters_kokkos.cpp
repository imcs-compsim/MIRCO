#include "mirco_inputparameters_kokkos.h"

#include <map>

#include "mirco_topology_kokkos.h"

namespace
{
  // Shape factors (See section 3.3 of https://doi.org/10.1007/s00466-019-01791-3)
  // These are the shape factors to calculate the elastic compliance correction of the micro-scale
  // contact constitutive law for various resolutions.
  // NOTE: Currently MIRCO works for resoulution of 1 to 8. The following map store the shape
  // factors for resolution of 1 to 8.

  // The following pressure based constants are calculated by solving a flat indentor problem using
  // the pressure based Green function described in Pohrt and Li (2014).
  // http://dx.doi.org/10.1134/s1029959914040109
  const std::map<int, double> shape_factors_pressure{{1, 0.961389237917602}, {2, 0.924715342432435},
      {3, 0.899837531880697}, {4, 0.884976751041942}, {5, 0.876753783192863},
      {6, 0.872397956576882}, {7, 0.8701463093314326}, {8, 0.8689982669426167}};

  // The following force based constants are taken from Table 1 of Bonari et al. (2020).
  // https://doi.org/10.1007/s00466-019-01791-3
  const std::map<int, double> shape_factors_force{{1, 0.778958541513360}, {2, 0.805513388666376},
      {3, 0.826126871395416}, {4, 0.841369158110513}, {5, 0.851733020725652},
      {6, 0.858342234203154}, {7, 0.862368243479785}, {8, 0.864741597831785}};

  double InterpolatedShapeFactor(const std::map<int, double>& shapeFactors, int N)
  {
    const double resolution = log2(N - 1);
    const int resFloor = static_cast<int>(std::floor(resolution));
    const double sfFloor = shapeFactors.at(resFloor);
    if (((N - 1) & (N - 2)) == 0) return sfFloor;
    return sfFloor + (resolution - resFloor) * (shapeFactors.at(resFloor + 1) - sfFloor);
  }
}  // namespace

namespace MIRCO
{
  InputParameters::InputParameters(double E1, double E2, double nu1, double nu2, double Tolerance,
      double Delta, double LateralLength, int Resolution, double InitialTopologyStdDeviation,
      double Hurst, bool RandomSeedFlag, int RandomGeneratorSeed, int MaxIteration,
      bool WarmStartingFlag, bool PressureGreenFunFlag)
      : tolerance(Tolerance),
        delta(Delta),
        lateral_length(LateralLength),
        max_iteration(MaxIteration),
        warm_starting_flag(WarmStartingFlag),
        pressure_green_funct_flag(PressureGreenFunFlag)
  {
    auto topology_h = CreateRmgSurface(
        Resolution, InitialTopologyStdDeviation, Hurst, RandomSeedFlag, RandomGeneratorSeed);
    topology_d = Kokkos::create_mirror_view_and_copy(ExecSpace_Default_t(), topology_h);

    // resolution is available; no interpolation needed
    if (PressureGreenFunFlag)
    {
      shape_factor = shape_factors_pressure.at(Resolution);
    }
    else
    {
      shape_factor = shape_factors_force.at(Resolution);
    }

    composite_youngs = 1.0 / ((1 - nu1 * nu1) / E1 + (1 - nu2 * nu2) / E2);
    elastic_compliance_correction = LateralLength * composite_youngs / shape_factor;
    N = (1 << Resolution) + 1;
    grid_size = LateralLength / N;
  }

  InputParameters::InputParameters(double E1, double E2, double nu1, double nu2, double Tolerance,
      double Delta, double LateralLength, const std::string& TopologyFilePath, int MaxIteration,
      bool WarmStartingFlag, bool PressureGreenFunFlag)
      : tolerance(Tolerance),
        delta(Delta),
        lateral_length(LateralLength),
        max_iteration(MaxIteration),
        warm_starting_flag(WarmStartingFlag),
        pressure_green_funct_flag(PressureGreenFunFlag)
  {
    auto topology_h = CreateSurfaceFromFile(TopologyFilePath, N);
    topology_d = Kokkos::create_mirror_view_and_copy(ExecSpace_Default_t(), topology_h);

    // interpolation needed
    if (PressureGreenFunFlag)
    {
      shape_factor = InterpolatedShapeFactor(shape_factors_pressure, N);
    }
    else
    {
      shape_factor = InterpolatedShapeFactor(shape_factors_force, N);
    }

    composite_youngs = 1.0 / ((1 - nu1 * nu1) / E1 + (1 - nu2 * nu2) / E2);
    elastic_compliance_correction = LateralLength * composite_youngs / shape_factor;
    grid_size = LateralLength / N;
  }

}  // namespace MIRCO
