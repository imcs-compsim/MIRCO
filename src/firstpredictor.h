#ifndef SRC_FIRSTPREDICTOR_H_
#define SRC_FIRSTPREDICTOR_H_

#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>

void FirstPredictor(int &n0, std::vector<double> &xv0, std::vector<double> &yv0,
    std::vector<double> &b0, double zmax, double Delta, double w_el, std::vector<double> x,
    Epetra_SerialDenseMatrix topology);

#endif  // SRC_FIRSTPREDICTOR_H_