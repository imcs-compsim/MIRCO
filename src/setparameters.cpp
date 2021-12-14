
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <jsoncpp/json/json.h>
using namespace std;

#include "filesystem_utils.h"
#include "setparameters.h"

void SetParameters(double& E1, double& E2,
                   double& lato, double& nu1,
                   double& nu2, double& G1, double& G2, double& E,
                   double& alpha,
                   double& k_el, double& delta, double& nnodi, double& errf,
                   double& to1, double& Delta, string& zfilePath, int& n, string& jsonFileName, bool& rmg_flag, double& Hurst, bool& rand_seed_flag) {
  
  
  Json::Value parameterlist;   // will contain the root value after parsing.
  ifstream stream(jsonFileName, std::ifstream::binary);
  stream >> parameterlist; 

  rmg_flag = parameterlist["rmg_flag"].asBool();
  rand_seed_flag = parameterlist["rand_seed_flag"].asBool();
  zfilePath = parameterlist["z_file_path"].asString();

  ChangeRelativePath(zfilePath, jsonFileName);

  E1 = parameterlist["parameters"]["material_parameters"]["E1"].asDouble();
  E2 = parameterlist["parameters"]["material_parameters"]["E2"].asDouble();
  nu1 = parameterlist["parameters"]["material_parameters"]["nu1"].asDouble();
  nu2 = parameterlist["parameters"]["material_parameters"]["nu2"].asDouble();
  G1 = E1 / (2 * (1 + nu1));
  G2 = E2 / (2 * (1 + nu2));
  E = 1 / ((1 - pow(nu1, 2)) / E1 + (1 - pow(nu2, 2) / E2));
  vector<double> alpha_con{0.778958541513360, 0.805513388666376,
                           0.826126871395416, 0.841369158110513,
                           0.851733020725652, 0.858342234203154,
                           0.862368243479785, 0.864741597831785};
  n = parameterlist["parameters"]["geometrical_parameters"]["n"].asInt();
  Hurst = parameterlist["parameters"]["geometrical_parameters"]["H"].asDouble(); // Hurst component
  alpha = alpha_con[n - 1];
  lato = parameterlist["parameters"]["geometrical_parameters"]["lato"].asDouble();  // Lateral side of the surface [micrometers]
  k_el = lato * E / alpha; 
  delta = lato / (pow(2, n) + 1);
  nnodi = pow(pow(2, n + 1), 2);
  errf = parameterlist["parameters"]["geometrical_parameters"]["errf"].asDouble();
  to1 = parameterlist["parameters"]["geometrical_parameters"]["tol"].asDouble();
  Delta = parameterlist["parameters"]["geometrical_parameters"]["Delta"].asDouble();
}