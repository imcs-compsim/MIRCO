
#include <string>
#include <vector>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "filesystem_utils.h"
#include "setparameters.h"

void SetParameters(double& E1, double& E2, double& lato, double& nu1, double& nu2, double& G1,
    double& G2, double& E, double& alpha, double& k_el, double& delta, double& nnodi, double& errf,
    double& tol, double& Delta, std::string& zfilePath, int& n, const std::string& inputFileName,
    bool& rmg_flag, double& Hurst, bool& rand_seed_flag, int& rmg_seed, bool& flagwarm)
{
  Teuchos::RCP<Teuchos::ParameterList> parameterList = Teuchos::rcp(new Teuchos::ParameterList());
  Teuchos::updateParametersFromXmlFile(inputFileName, parameterList.ptr());

  flagwarm = parameterList->get<bool>("flagwarm");
  rmg_flag = parameterList->get<bool>("rmg_flag");
  rand_seed_flag = parameterList->get<bool>("rand_seed_flag");
  rmg_seed = parameterList->get<int>("rmg_seed");
  zfilePath = parameterList->get<std::string>("z_file_path");

  UTILS::ChangeRelativePath(zfilePath, inputFileName);

  Teuchos::ParameterList& matParams =
      parameterList->sublist("parameters").sublist("material_parameters");
  E1 = matParams.get<double>("E1");
  E2 = matParams.get<double>("E2");
  nu1 = matParams.get<double>("nu1");
  nu2 = matParams.get<double>("nu2");
  G1 = E1 / (2 * (1 + nu1));
  G2 = E2 / (2 * (1 + nu2));
  E = 1 / ((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2) / E2));
  std::vector<double> alpha_con{0.778958541513360, 0.805513388666376, 0.826126871395416,
      0.841369158110513, 0.851733020725652, 0.858342234203154, 0.862368243479785,
      0.864741597831785};

  Teuchos::ParameterList& geoParams =
      parameterList->sublist("parameters").sublist("geometrical_parameters");
  n = geoParams.get<int>("n");
  Hurst = geoParams.get<double>("H");    // Hurst exponent
  lato = geoParams.get<double>("lato");  // Lateral side of the surface [micrometers]
  errf = geoParams.get<double>("errf");
  tol = geoParams.get<double>("tol");
  Delta = geoParams.get<double>("Delta");

  alpha = alpha_con[n - 1];
  k_el = lato * E / alpha;
  delta = lato / (pow(2, n) + 1);
  nnodi = pow(pow(2, n + 1), 2);
}
