#ifndef SRC_SECONDPREDICTOR_H_
#define SRC_SECONDPREDICTOR_H_

#include <vector>

void SecondPredictor(bool flagwarm, int k, int n0, int nf2, std::vector<double> xv0,
    std::vector<double> yv0, std::vector<double> pf, std::vector<double> &x0,
    std::vector<double> &b0, std::vector<double> xvf, std::vector<double> yvf);

#endif  // SRC_SECONDPREDICTOR_H_