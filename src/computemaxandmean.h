#ifndef SRC_COMPUTEMAXANDMEAN_H_
#define SRC_COMPUTEMAXANDMEAN_H_

#include <Epetra_SerialSymDenseMatrix.h>

void ComputeMaxAndMean(Epetra_SerialDenseMatrix topology, double &zmax, double &zmean);

#endif  // SRC_COMPUTEMAXANDMEAN_H_