
#include <string>
#include <vector>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "filesystem_utils.h"
#include "setparameters.h"

void MIRCO::SetParameters(double& E1, double& E2, double& lato, double& nu1, double& nu2,
    double& G1, double& G2, double& E, double& alpha, double& k_el, double& delta, double& nnodi,
    double& errf, double& tol, double& Delta, std::string& zfilePath, int& resolution,
    const std::string& inputFileName, bool& rmg_flag, double& Hurst, bool& rand_seed_flag,
    int& rmg_seed, bool& flagwarm, int& max_iter)
{
  // The aim of this fuction is to read the input file containing the simulation specific, material
  // and geometrical parameters.

  Teuchos::RCP<Teuchos::ParameterList> parameterList = Teuchos::rcp(new Teuchos::ParameterList());
  Teuchos::updateParametersFromXmlFile(inputFileName, parameterList.ptr());

  // Setting up the simulation specific parameters.

  // set flagwarm --> true for using the warm starter. It predicts the nodes coming into contact in
  // the next iteration and hence speeds up the computation.
  flagwarm = parameterList->get<bool>("flagwarm");
  // set rmg_flag --> true to use the random mid-point generator to generate a topology.
  // set rmg_flag --> flase to read the topology from an input file.
  rmg_flag = parameterList->get<bool>("rmg_flag");
  // set rand_sed_flag --> true to fix the seed to generate psuedo random topology to reproduce
  // results.
  rand_seed_flag = parameterList->get<bool>("rand_seed_flag");
  // rmg_seed --> set the value of seed for the random mid-point generator
  rmg_seed = parameterList->get<int>("rmg_seed");
  // zfilePath --> name of the input file containing the topology.
  zfilePath = parameterList->get<std::string>("z_file_path");
  // max_iter --> maximum number of iterations for the force to converge.
  max_iter = parameterList->get<int>("max_iter");

  // following function generates the actual path of the topology file.
  UTILS::ChangeRelativePath(zfilePath, inputFileName);

  // Setting up the material parameters.
  Teuchos::ParameterList& matParams =
      parameterList->sublist("parameters").sublist("material_parameters");

  // E1 and E2 are the Young's moduli of bodies 1 and 2 respectively.
  E1 = matParams.get<double>("E1");
  E2 = matParams.get<double>("E2");
  // nu1 and nu2 are the Poisson's ratios.
  nu1 = matParams.get<double>("nu1");
  nu2 = matParams.get<double>("nu2");
  // G1 and G2 are respective Shear moduli.
  G1 = E1 / (2 * (1 + nu1));
  G2 = E2 / (2 * (1 + nu2));
  // E is the composite Young's modulus.
  E = pow(((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2)) / E2), -1);
  // alpha_con is a vector containiing the shape factor for various resolution parameters.
  std::vector<double> alpha_con{0.778958541513360, 0.805513388666376, 0.826126871395416,
      0.841369158110513, 0.851733020725652, 0.858342234203154, 0.862368243479785,
      0.864741597831785};

  // Setting up the geometrical parameters.
  Teuchos::ParameterList& geoParams =
      parameterList->sublist("parameters").sublist("geometrical_parameters");
  // n --> resolution parameter
  resolution = geoParams.get<int>("n");
  // Hurst --> Hurst Exponent (Used in random mid-point generator)
  Hurst = geoParams.get<double>("H");
  // lato --> Lateral side of the surface [micrometers]
  lato = geoParams.get<double>("lato");
  // errf --> initialing the error in the force vector.
  errf = geoParams.get<double>("errf");
  // tol --> tolerance for the convergence of force.
  tol = geoParams.get<double>("tol");
  // Delta --> far-field displacement.
  Delta = geoParams.get<double>("Delta");

  // alpha --> shape factor for this problem
  alpha = alpha_con[resolution - 1];
  // k_el --> Elastic compliance correction
  k_el = lato * E / alpha;
  // delta --> grid size
  delta = lato / (pow(2, resolution) + 1);
  // nnodi --> total number of nodes
  nnodi = pow(pow(2, resolution + 1), 2);
}
