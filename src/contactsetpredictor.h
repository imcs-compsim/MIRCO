#ifndef SRC_CONTACTSETPREDICTOR_H_
#define SRC_CONTACTSETPREDICTOR_H_

#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>

void ContactSetPredictor(int &n0, std::vector<double> &xv0, std::vector<double> &yv0,
    std::vector<double> &b0, double zmax, double Delta, double w_el, std::vector<double> x,
    Epetra_SerialDenseMatrix topology);

#endif  // SRC_CONTACTSETPREDICTOR_H_