#ifndef SRC_INITIALGUESSPREDICTOR_H_
#define SRC_INITIALGUESSPREDICTOR_H_

#include <vector>

void InitialGuessPredictor(bool flagwarm, int k, int n0, int nf2, std::vector<double> xv0,
    std::vector<double> yv0, std::vector<double> pf, std::vector<double> &x0,
    std::vector<double> &b0, std::vector<double> xvf, std::vector<double> yvf);

#endif  // SRC_INITIALGUESSPREDICTOR_H_