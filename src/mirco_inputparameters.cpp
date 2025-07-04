#include "mirco_inputparameters.h"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <map>
#include <string>
#include <vector>

#include "mirco_filesystem_utils.h"
#include "mirco_topology.h"
#include "mirco_topologyutilities.h"

inline double Resolution_wrtN(int N) { return log2(N - 1); }
inline int N_wrtResolution(int resolution)
{
  return (1 << resolution) + 1;  // 2^resolution + 1
}

// basic linear interpolation for now
double InterpolatedShapeFactor_wrtN(std::map<int, double>& shapeFactors, int N)
{
  double resolution = Resolution_wrtN(N);
  int floor = (int)std::floor(resolution);
  double sfFloor = shapeFactors.at(floor);
  return sfFloor + (resolution - floor) * (shapeFactors.at(floor + 1) - sfFloor);
}

MIRCO::InputParameters::InputParameters(const std::string& inputFileName)
{
  Teuchos::RCP<Teuchos::ParameterList> parameter_list = Teuchos::rcp(new Teuchos::ParameterList());
  Teuchos::updateParametersFromXmlFile(inputFileName, parameter_list.ptr());

  // Setting up the simulation specific parameters.
  warm_starting_flag_ = parameter_list->get<bool>("WarmStartingFlag");


  max_iteration_ = parameter_list->get<int>("MaxIteration");

  // Setting up the geometrical parameters.
  Teuchos::ParameterList& geo_params =
      parameter_list->sublist("parameters").sublist("geometrical_parameters");
  lateral_length_ = geo_params.get<double>("LateralLength");
  tolerance_ = geo_params.get<double>("Tolerance");
  delta_ = geo_params.get<double>("Delta");

  // Shape factors (See section 3.3 of https://doi.org/10.1007/s00466-019-01791-3)
  // These are the shape factors to calculate the elastic compliance correction of the micro-scale
  // contact constitutive law for various resolutions.
  // NOTE: Currently MIRCO works for resouluion of 1 to 8. The following map store the shape
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

  // Set the surface generator based on RandomTopologyFlag
  bool random_topology_flag = parameter_list->get<bool>("RandomTopologyFlag");
  if (random_topology_flag)
  {
    auto resolution = geo_params.get<int>("Resolution");

    auto surfacegenerator = Teuchos::rcp(
        new MIRCO::Rmg(resolution, geo_params.get<double>("InitialTopologyStdDeviation"),
            geo_params.get<double>("HurstExponent"), parameter_list->get<bool>("RandomSeedFlag"),
            parameter_list->get<int>("RandomGeneratorSeed")));



    // resolution is available; no interpolation needed
    if (parameter_list->get<bool>("PressureGreenFunFlag"))
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
    auto topology_file_path = parameter_list->get<std::string>("TopologyFilePath");
    // following function generates the actual path of the topology file.
    MIRCO::UTILS::ChangeRelativePath(topology_file_path, inputFileName);

    auto surface_generator = Teuchos::rcp(new MIRCO::ReadFile(topology_file_path));

    // interpolation needed
    if (parameter_list->get<bool>("PressureGreenFunFlag"))
    {
      shape_factor_ = InterpolatedShapeFactors_wrtN(shape_factors_pressure, )
    }
    else
    {
      shape_factor_ = shape_factors_force.at(resolution);
    }
  }



  // Setting up the material parameters.
  Teuchos::ParameterList& matParams =
      parameter_list->sublist("parameters").sublist("material_parameters");

  double E1 = matParams.get<double>("E1");
  double E2 = matParams.get<double>("E2");
  double nu1 = matParams.get<double>("nu1");
  double nu2 = matParams.get<double>("nu2");

  // Composite Young's modulus.
  composite_youngs_ = pow(((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2)) / E2), -1);



  elastic_compliance_correction_ = lateral_length_ * composite_youngs_ / shape_factor_;
  grid_size_ = lateral_length_ / (pow(2, Resolution) + 1);
}



/*
Evaluate(pressure, Delta, LateralLength, GridSize, Tolerance, MaxIteration,
      CompositeYoungs, WarmStartingFlag, ElasticComplianceCorrection, topology, max_and_mean.max_,
      meshgrid, PressureGreenFunFlag
      */
