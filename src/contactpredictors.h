#ifndef SRC_CONTACTPREDICTORS_H_
#define SRC_CONTACTPREDICTORS_H_

#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>

void ContactSetPredictor(int &n0, std::vector<double> &xv0, std::vector<double> &yv0,
    std::vector<double> &b0, double zmax, double Delta, double w_el, std::vector<double> x,
    Epetra_SerialDenseMatrix topology);

void InitialGuessPredictor(bool flagwarm, int k, int n0, int nf2, std::vector<double> xv0,
    std::vector<double> yv0, std::vector<double> pf, std::vector<double> &x0,
    std::vector<double> &b0, std::vector<double> xvf, std::vector<double> yvf);

#endif  // SRC_CONTACTPREDICTORS_H_