
#include <string>
#include <vector>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "mirco_filesystem_utils.h"
#include "mirco_setparameters.h"

void MIRCO::SetParameters(double& E1, double& E2, double& LateralLength, double& nu1, double& nu2,
    double& CompositeYoungs, double& alpha, double& ElasticComplianceCorrection, double& GridSize,
    double& Tolerance, double& Delta, std::string& TopologyFilePath, int& Resolution,
    double& InitialStdDeviation, const std::string& inputFileName, bool& RandomTopologyFlag,
    double& Hurst, bool& RandomSeedFlag, int& RandomGeneratorSeed, bool& WarmStartingFlag,
    int& MaxIteration)
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

  // Correction factor vector
  std::vector<double> alpha_con{0.778958541513360, 0.805513388666376, 0.826126871395416,
      0.841369158110513, 0.851733020725652, 0.858342234203154, 0.862368243479785,
      0.864741597831785};

  // Setting up the geometrical parameters.
  Teuchos::ParameterList& geoParams =
      parameterList->sublist("parameters").sublist("geometrical_parameters");

  Resolution = geoParams.get<int>("Resolution");
  Hurst = geoParams.get<double>("HurstExponent");
  LateralLength = geoParams.get<double>("LateralLength");
  InitialStdDeviation = geoParams.get<double>("InitialStdDeviation");
  Tolerance = geoParams.get<double>("Tolerance");
  Delta = geoParams.get<double>("Delta");
  alpha = alpha_con[Resolution - 1];
  ElasticComplianceCorrection = LateralLength * CompositeYoungs / alpha;
  GridSize = LateralLength / (pow(2, Resolution) + 1);
}
