
#include "mirco_setparameters.h"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <map>
#include <string>
#include <vector>

#include "mirco_filesystem_utils.h"

void MIRCO::SetParameters(double& E1, double& E2, double& LateralLength, double& nu1, double& nu2,
    double& CompositeYoungs, double& ShapeFactor, double& ElasticComplianceCorrection,
    double& GridSize, double& Tolerance, double& Delta, std::string& TopologyFilePath,
    int& Resolution, double& InitialTopologyStdDeviation, const std::string& inputFileName,
    bool& RandomTopologyFlag, double& Hurst, bool& RandomSeedFlag, int& RandomGeneratorSeed,
    bool& WarmStartingFlag, int& MaxIteration, bool& PressureGreenFunFlag)
{
  Teuchos::RCP<Teuchos::ParameterList> parameterList = Teuchos::rcp(new Teuchos::ParameterList());
  Teuchos::updateParametersFromXmlFile(inputFileName, parameterList.ptr());

  // Setting up the simulation specific parameters.
  WarmStartingFlag = parameterList->get<bool>("WarmStartingFlag");
  RandomTopologyFlag = parameterList->get<bool>("RandomTopologyFlag");
  RandomSeedFlag = parameterList->get<bool>("RandomSeedFlag");
  RandomGeneratorSeed = parameterList->get<int>("RandomGeneratorSeed");
  TopologyFilePath = parameterList->get<std::string>("TopologyFilePath");
  MaxIteration = parameterList->get<int>("MaxIteration");
  PressureGreenFunFlag = parameterList->get<bool>("PressureGreenFunFlag");

  // following function generates the actual path of the topology file.
  MIRCO::UTILS::ChangeRelativePath(TopologyFilePath, inputFileName);

  // Setting up the material parameters.
  Teuchos::ParameterList& matParams =
      parameterList->sublist("parameters").sublist("material_parameters");

  E1 = matParams.get<double>("E1");
  E2 = matParams.get<double>("E2");
  nu1 = matParams.get<double>("nu1");
  nu2 = matParams.get<double>("nu2");

  // Composite Young's modulus.
  CompositeYoungs = pow(((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2)) / E2), -1);

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


  // Setting up the geometrical parameters.
  Teuchos::ParameterList& geoParams =
      parameterList->sublist("parameters").sublist("geometrical_parameters");

  Resolution = geoParams.get<int>("Resolution");
  Hurst = geoParams.get<double>("HurstExponent");
  LateralLength = geoParams.get<double>("LateralLength");
  InitialTopologyStdDeviation = geoParams.get<double>("InitialTopologyStdDeviation");
  Tolerance = geoParams.get<double>("Tolerance");
  Delta = geoParams.get<double>("Delta");

  if (PressureGreenFunFlag)
  {
    ShapeFactor = shape_factors_pressure.at(Resolution);
  }
  else
  {
    ShapeFactor = shape_factors_force.at(Resolution);
  }

  ElasticComplianceCorrection = LateralLength * CompositeYoungs / ShapeFactor;
  GridSize = LateralLength / (pow(2, Resolution) + 1);
}
