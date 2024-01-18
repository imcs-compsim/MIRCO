
#include "mirco_setparameters.h"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <string>
#include <vector>

#include "mirco_filesystem_utils.h"

void MIRCO::SetParameters(double& E1, double& E2, double& LateralLength, double& nu1, double& nu2,
    double& CompositeYoungs, double& CompositePoissonsRatio, double& alpha,
    double& ElasticComplianceCorrection, double& GridSize, double& Tolerance, double& Delta,
    std::string& TopologyFilePath, int& Resolution, double& InitialTopologyStdDeviation,
    const std::string& inputFileName, bool& RandomTopologyFlag, double& Hurst, bool& RandomSeedFlag,
    int& RandomGeneratorSeed, bool& WarmStartingFlag, int& MaxIteration, bool& PressureGreenFunFlag)
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
  UTILS::ChangeRelativePath(TopologyFilePath, inputFileName);

  // Setting up the material parameters.
  Teuchos::ParameterList& matParams =
      parameterList->sublist("parameters").sublist("material_parameters");

  E1 = matParams.get<double>("E1");
  E2 = matParams.get<double>("E2");
  nu1 = matParams.get<double>("nu1");
  nu2 = matParams.get<double>("nu2");

  // Composite Young's modulus.
  CompositeYoungs = pow(((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2)) / E2), -1);

  // Composite Shear modulus
  double G1 = E1 / (2 * (1 + nu1));
  double G2 = E2 / (2 * (1 + nu2));
  double CompositeShear = pow(((2 - nu1) / (4 * G1) + (2 - nu2) / (4 * G2)), -1);

  // Composite Poisson's ratio
  CompositePoissonsRatio = CompositeYoungs / (2 * CompositeShear) - 1;

  // Correction factor vectors
  // These are the correction factors to calculate the elastic compliance of the micro-scale contact
  // constitutive law for various resolutions.
  // The following pressure based constants are calculated by solving a flat indentor problem using
  // the pressure based Green function described in Pohrt and Li (2014).
  // http://dx.doi.org/10.1134/s1029959914040109
  std::vector<double> alpha_con_pressure{0.961389237917602, 0.924715342432435, 0.899837531880697,
      0.884976751041942, 0.876753783192863, 0.872397956576882, 0.871958228537090,
      0.882669916668780};
  // The following force based constants are taken from Table 1 of Bonari et al. (2020).
  // https://doi.org/10.1007/s00466-019-01791-3
  std::vector<double> alpha_con_force{0.778958541513360, 0.805513388666376, 0.826126871395416,
      0.841369158110513, 0.851733020725652, 0.858342234203154, 0.862368243479785,
      0.864741597831785};

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
    alpha = alpha_con_pressure[Resolution - 1];
  }
  else
  {
    alpha = alpha_con_force[Resolution - 1];
  }

  ElasticComplianceCorrection = LateralLength * CompositeYoungs / alpha;
  GridSize = LateralLength / (pow(2, Resolution) + 1);
}
