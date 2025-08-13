#include "mirco_inputparameters_kokkos.h"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <map>
#include <string>

#include "mirco_filesystem_utils.h"
#include "mirco_topology_kokkos.h"

// Shape factors (See section 3.3 of https://doi.org/10.1007/s00466-019-01791-3)
// These are the shape factors to calculate the elastic compliance correction of the micro-scale
// contact constitutive law for various resolutions.
// NOTE: Currently MIRCO works for resoulution of 1 to 8. The following map store the shape
// factors for resolution of 1 to 8.

// The following pressure based constants are calculated by solving a flat indentor problem using
// the pressure based Green function described in Pohrt and Li (2014).
// http://dx.doi.org/10.1134/s1029959914040109
const std::map<int, double> shape_factors_pressure{{1, 0.961389237917602}, {2, 0.924715342432435},
    {3, 0.899837531880697}, {4, 0.884976751041942}, {5, 0.876753783192863}, {6, 0.872397956576882},
    {7, 0.8701463093314326}, {8, 0.8689982669426167}};

// The following force based constants are taken from Table 1 of Bonari et al. (2020).
// https://doi.org/10.1007/s00466-019-01791-3
const std::map<int, double> shape_factors_force{{1, 0.778958541513360}, {2, 0.805513388666376},
    {3, 0.826126871395416}, {4, 0.841369158110513}, {5, 0.851733020725652}, {6, 0.858342234203154},
    {7, 0.862368243479785}, {8, 0.864741597831785}};

inline int getNFromResolution(int resolution) { return (1 << resolution) + 1; }
double InterpolatedShapeFactor(const std::map<int, double>& shapeFactors, int N)
{
  const double resolution = log2(N - 1);
  const int resFloor = static_cast<int>(std::floor(resolution));
  const double sfFloor = shapeFactors.at(resFloor);
  if (((N - 1) & (N - 2)) == 0) return sfFloor;
  return sfFloor + (resolution - resFloor) * (shapeFactors.at(resFloor + 1) - sfFloor);
}

MIRCO::InputParameters::InputParameters(const std::string& inputFileName)
{
  Teuchos::ParameterList parameter_list;
  Teuchos::updateParametersFromXmlFile(inputFileName, Teuchos::ptrFromRef(parameter_list));

  Teuchos::ParameterList& geo_params =
      parameter_list.sublist("parameters").sublist("geometrical_parameters");
  Teuchos::ParameterList& matParams =
      parameter_list.sublist("parameters").sublist("material_parameters");

  // Set the surface generator based on RandomTopologyFlag
  if (parameter_list.get<bool>("RandomTopologyFlag"))
  {
    *this = InputParameters(matParams.get<double>("E1"), matParams.get<double>("E2"),
        matParams.get<double>("nu1"), matParams.get<double>("nu2"),
        geo_params.get<double>("Tolerance"), geo_params.get<double>("Delta"),
        geo_params.get<double>("LateralLength"), geo_params.get<int>("Resolution"),
        geo_params.get<double>("InitialTopologyStdDeviation"),
        geo_params.get<double>("HurstExponent"), parameter_list.get<bool>("RandomSeedFlag"),
        parameter_list.get<int>("RandomGeneratorSeed"), parameter_list.get<int>("MaxIteration"),
        parameter_list.get<bool>("WarmStartingFlag"),
        parameter_list.get<bool>("PressureGreenFunFlag"));
  }
  else
  {
    auto topology_file_path = parameter_list.get<std::string>("TopologyFilePath");
    // The following function generates the actual path of the topology file
    MIRCO::UTILS::ChangeRelativePath(topology_file_path, inputFileName);

    *this = InputParameters(matParams.get<double>("E1"), matParams.get<double>("E2"),
        matParams.get<double>("nu1"), matParams.get<double>("nu2"),
        geo_params.get<double>("Tolerance"), geo_params.get<double>("Delta"),
        geo_params.get<double>("LateralLength"), topology_file_path,
        parameter_list.get<int>("MaxIteration"), parameter_list.get<bool>("WarmStartingFlag"),
        parameter_list.get<bool>("PressureGreenFunFlag"));
  }
}

MIRCO::InputParameters::InputParameters(double E1, double E2, double nu1, double nu2,
    double Tolerance, double Delta, double LateralLength, int Resolution,
    double InitialTopologyStdDeviation, double Hurst, bool RandomSeedFlag, int RandomGeneratorSeed,
    int MaxIteration, bool WarmStartingFlag, bool PressureGreenFunFlag)
    : tolerance(Tolerance),
      delta(Delta),
      lateral_length(LateralLength),
      max_iteration(MaxIteration),
      warm_starting_flag(WarmStartingFlag),
      pressure_green_funct_flag(PressureGreenFunFlag)
{
  auto topology_h = MIRCO::CreateRmgSurface(
      Resolution, InitialTopologyStdDeviation, Hurst, RandomSeedFlag, RandomGeneratorSeed);
  topology_d = Kokkos::create_mirror_view_and_deep_copy(topology_h);

  // resolution is available; no interpolation needed
  if (PressureGreenFunFlag)
  {
    shape_factor = shape_factors_pressure.at(Resolution);
  }
  else
  {
    shape_factor = shape_factors_force.at(Resolution);
  }

  composite_youngs = pow(((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2)) / E2), -1);
  elastic_compliance_correction = LateralLength * composite_youngs / shape_factor;
  N = NFromResolution(Resolution);
  grid_size = LateralLength / N;
}

MIRCO::InputParameters::InputParameters(double E1, double E2, double nu1, double nu2,
    double Tolerance, double Delta, double LateralLength, const std::string& TopologyFilePath,
    int MaxIteration, bool WarmStartingFlag, bool PressureGreenFunFlag)
    : tolerance(Tolerance),
      delta(Delta),
      lateral_length(LateralLength),
      max_iteration(MaxIteration),
      warm_starting_flag(WarmStartingFlag),
      pressure_green_funct_flag(PressureGreenFunFlag)
{
  auto topology_h = std::make_shared<Teuchos::SerialDenseMatrix<int, double>>(
      MIRCO::CreateSurfaceFromFile(TopologyFilePath, N));
  topology_d = Kokkos::create_mirror_view_and_deep_copy(topology_h);

  // interpolation needed
  if (PressureGreenFunFlag)
  {
    shape_factor = InterpolatedShapeFactor(shape_factors_pressure, N);
  }
  else
  {
    shape_factor = InterpolatedShapeFactor(shape_factors_force, N);
  }

  composite_youngs = pow(((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2)) / E2), -1);
  elastic_compliance_correction = LateralLength * composite_youngs / shape_factor;
  grid_size = LateralLength / N;
}
