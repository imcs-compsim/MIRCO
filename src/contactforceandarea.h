#ifndef SRC_CONTACTFORCEANDAREA_H_
#define SRC_CONTACTFORCEANDAREA_H_

#include <cmath>
#include <vector>

void ContactForceAndArea(std::vector<double> &force0, std::vector<double> &area0, int &iter,
    double &w_el, double nf, std::vector<double> pf, int k, double delta, double lato, double k_el);

#endif  // SRC_CONTACTFORCEANDAREA_H_