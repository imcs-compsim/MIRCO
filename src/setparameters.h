#ifndef SRC_SETPARAMETERS_H_
#define SRC_SETPARAMETERS_H_

#include <string>

void SetParameters(double& E1, double& E2, double& lato, double& nu1, double& nu2, double& G1,
    double& G2, double& E, double& alpha, double& k_el, double& delta, double& nnodi, double& errf,
    double& tol, double& Delta, std::string& zfilePath, int& n, const std::string& inputFileName,
    bool& rmg_flag, double& Hurst, bool& rand_seed_flag, int& rmg_seed, bool& flagwarm,
    int& max_iter);


#endif  // SRC_SETPARAMETERS_H_
