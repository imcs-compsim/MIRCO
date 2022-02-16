#ifndef SRC_SURFMAXANDMEAN_H_
#define SRC_SURFMAXANDMEAN_H_

#include <Epetra_SerialSymDenseMatrix.h>

void SurfMaxAndMean(Epetra_SerialDenseMatrix topology, double &zmax, double &zmean);

#endif  // SRC_SURFMAXANDMEAN_H_