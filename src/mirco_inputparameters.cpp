#include "mirco_inputparameters.h"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <map>
#include <string>

#include "mirco_filesystem_utils.h"
#include "mirco_topology.h"

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

inline int N_wrtResolution(int resolution)
{
  // TODO: Add checks
  return (1 << resolution) + 1;  // 2^resolution + 1
}
// basic linear interpolation for now
double InterpolatedShapeFactor_wrtN(const std::map<int, double>& shapeFactors, int N)
{
  double resolution = log2(N - 1);
  int floor = (int)std::floor(resolution);
  double sfFloor = shapeFactors.at(floor);
  return sfFloor + (resolution - floor) * (shapeFactors.at(floor + 1) - sfFloor);
}

MIRCO::InputParameters::InputParameters(const std::string& inputFileName)
{
  Teuchos::ParameterList parameter_list;
  Teuchos::updateParametersFromXmlFile(inputFileName, Teuchos::ptrFromRef(parameter_list));

  // Setting up the simulation specific parameters.
  warm_starting_flag_ = parameter_list.get<bool>("WarmStartingFlag");
  max_iteration_ = parameter_list.get<int>("MaxIteration");

  // Setting up the geometrical parameters.
  Teuchos::ParameterList& geo_params =
      parameter_list.sublist("parameters").sublist("geometrical_parameters");
  lateral_length_ = geo_params.get<double>("LateralLength");
  tolerance_ = geo_params.get<double>("Tolerance");
  delta_ = geo_params.get<double>("Delta");

  int N;
  // Set the surface generator based on RandomTopologyFlag
  if (parameter_list.get<bool>("RandomTopologyFlag"))
  {
    auto resolution = geo_params.get<int>("Resolution");

    N = N_wrtResolution(resolution);

    topology_ =
        MIRCO::CreateRmgSurface(resolution, geo_params.get<double>("InitialTopologyStdDeviation"),
            geo_params.get<double>("HurstExponent"), parameter_list.get<bool>("RandomSeedFlag"),
            parameter_list.get<int>("RandomGeneratorSeed"));



    // resolution is available; no interpolation needed
    if (parameter_list.get<bool>("PressureGreenFunFlag"))
    {
      shape_factor_ = shape_factors_pressure.at(resolution);
    }
    else
    {
      shape_factor_ = shape_factors_force.at(resolution);
    }
  }
  else
  {
    auto topology_file_path = parameter_list.get<std::string>("TopologyFilePath");
    // following function generates the actual path of the topology file.
    MIRCO::UTILS::ChangeRelativePath(topology_file_path, inputFileName);
    topology_ = MIRCO::CreateSurfaceFromFile(topology_file_path, N);

    // interpolation needed
    if (parameter_list.get<bool>("PressureGreenFunFlag"))
    {
      shape_factor_ = InterpolatedShapeFactor_wrtN(shape_factors_pressure, N);
    }
    else
    {
      shape_factor_ = InterpolatedShapeFactor_wrtN(shape_factors_force, N);
    }
  }



  // Setting up the material parameters.
  Teuchos::ParameterList& matParams =
      parameter_list.sublist("parameters").sublist("material_parameters");

  double E1 = matParams.get<double>("E1");
  double E2 = matParams.get<double>("E2");
  double nu1 = matParams.get<double>("nu1");
  double nu2 = matParams.get<double>("nu2");

  composite_youngs_ = pow(((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2)) / E2), -1);
  elastic_compliance_correction_ = lateral_length_ * composite_youngs_ / shape_factor_;
  grid_size_ = lateral_length_ / N;
}

MIRCO::InputParameters::InputParameters(double E1, double E2, double nu1, double nu2,
    double Tolerance, double Delta, int Resolution, double LateralLength,
    double InitialTopologyStdDeviation, double Hurst, bool RandomSeedFlag, int RandomGeneratorSeed,
    int MaxIteration, bool WarmStartingFlag, bool PressureGreenFunFlag)
    : tolerance_(Tolerance),
      delta_(Delta),
      lateral_length_(LateralLength),
      max_iteration_(MaxIteration),
      warm_starting_flag_(WarmStartingFlag),
      pressure_green_funct_flag_(PressureGreenFunFlag)
{
  topology_ = MIRCO::CreateRmgSurface(
      Resolution, InitialTopologyStdDeviation, Hurst, RandomSeedFlag, RandomGeneratorSeed);

  // resolution is available; no interpolation needed
  if (PressureGreenFunFlag)
  {
    shape_factor_ = shape_factors_pressure.at(Resolution);
  }
  else
  {
    shape_factor_ = shape_factors_force.at(Resolution);
  }

  composite_youngs_ = pow(((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2)) / E2), -1);
  elastic_compliance_correction_ = LateralLength * composite_youngs_ / shape_factor_;
  grid_size_ = LateralLength / N_wrtResolution(Resolution);
}

MIRCO::InputParameters::InputParameters(double E1, double E2, double nu1, double nu2,
    double Tolerance, double Delta, double LateralLength, const std::string& TopologyFilePath,
    int MaxIteration, bool WarmStartingFlag, bool PressureGreenFunFlag)
    : tolerance_(Tolerance),
      delta_(Delta),
      lateral_length_(LateralLength),
      max_iteration_(MaxIteration),
      warm_starting_flag_(WarmStartingFlag),
      pressure_green_funct_flag_(PressureGreenFunFlag)
{
  int N;
  topology_ = MIRCO::CreateSurfaceFromFile(TopologyFilePath, N);

  // interpolation needed
  if (PressureGreenFunFlag)
  {
    shape_factor_ = InterpolatedShapeFactor_wrtN(shape_factors_pressure, N);
  }
  else
  {
    shape_factor_ = InterpolatedShapeFactor_wrtN(shape_factors_force, N);
  }

  composite_youngs_ = pow(((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2)) / E2), -1);
  elastic_compliance_correction_ = LateralLength * composite_youngs_ / shape_factor_;
  grid_size_ = LateralLength / N;
}
