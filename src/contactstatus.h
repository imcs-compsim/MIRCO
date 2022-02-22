#ifndef SRC_CONTACTSTATUS_H_
#define SRC_CONTACTSTATUS_H_

#include <Epetra_SerialSymDenseMatrix.h>
#include <cmath>
#include <vector>

void ComputeContactNodes(std::vector<double> &xvf, std::vector<double> &yvf,
    std::vector<double> &pf, int &cont, double &nf, Epetra_SerialDenseMatrix y,
    std::vector<double> xv0, std::vector<double> yv0);

void ComputeContactForceAndArea(std::vector<double> &force0, std::vector<double> &area0, int &iter,
    double &w_el, double nf, std::vector<double> pf, int k, double delta, double lato, double k_el);

#endif  // SRC_CONTACTSTATUS_H_