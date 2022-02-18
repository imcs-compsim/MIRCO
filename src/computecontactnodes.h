#ifndef SRC_COMPUTECONTACTNODES_H_
#define SRC_COMPUTECONTACTNODES_H_

#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>

void ComputeContactNodes(std::vector<double> &xvf, std::vector<double> &yvf,
    std::vector<double> &pf, int &cont, double &nf, Epetra_SerialDenseMatrix y,
    std::vector<double> xv0, std::vector<double> yv0);

#endif  // SRC_COMPUTECONTACTNODES_H_