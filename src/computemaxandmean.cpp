#include "computemaxandmean.h"
#include <Epetra_SerialSymDenseMatrix.h>

void ComputeMaxAndMean(Epetra_SerialDenseMatrix topology, double &zmax, double &zmean)
{
#pragma omp parallel for schedule(guided, 16) reduction(+ : zmean) reduction(max : zmax)
  // Static and Guided seem even but Guided makes more sense
  for (int i = 0; i < topology.M(); i++)
  {
    for (int j = 0; j < topology.N(); j++)
    {
      zmean += topology(i, j);
      if (topology(i, j) > zmax)
      {
        zmax = topology(i, j);
      }
    }
  }

  zmean = zmean / (topology.N() * topology.M());
}