#include "warmstart.h"
#include <Epetra_SerialSymDenseMatrix.h>
#include <vector>

using namespace std;

Epetra_SerialDenseMatrix Warmstarter::Warmstart2(Epetra_SerialDenseMatrix xv0,
                                                 Epetra_SerialDenseMatrix yv0,
                                                 Epetra_SerialDenseMatrix &xvf,
                                                 Epetra_SerialDenseMatrix &yvf,
                                                 Epetra_SerialDenseMatrix &pf) {
  Epetra_SerialDenseMatrix x0;
  x0.Shape(xv0.N(), 1);
  Epetra_SerialDenseMatrix combinedMatrix;
  combinedMatrix.Shape(2, xv0.N());
  // matfin = [xv0, yv0]
#pragma omp parallel for schedule(static, 16) // Always same workload -> Static
  for (int i = 0; i < xv0.N(); i++) {
    combinedMatrix(0, i) = xv0(0, i);
  }

#pragma omp parallel for schedule(static, 16)
  for (int i = 0; i < yv0.N(); i++) {
    combinedMatrix(1, i) = yv0(0, i);
  }

  vector<int> index;
#pragma omp parallel for schedule(dynamic,                                     \
                                  16) // Workload can differ vastly -> Dynamic
  for (int i = 0; i < pf.N(); i++) {
    // ind=find(matfin(:,1)==xvf(i) & matfin(:,2)==yvf(i));
    for (int j = 0; j < xvf.N(); j++) {
      if ((combinedMatrix(j, 0) == xvf(i, 0)) &&
          (combinedMatrix(j, 1) == yvf(i, 0))) {
        index.push_back(j);
      }
    }

    // x0(ind,1)=pf(i);
#pragma omp parallel for schedule(static, 16) // Always same workload -> Static
    for (long unsigned int y = 0; y < index.size(); y++) {
      x0(y, 0) = pf(y, 0);
    }
    index.clear();
  }
  return x0;
}