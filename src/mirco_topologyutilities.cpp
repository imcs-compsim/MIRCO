#include "mirco_topologyutilities.h"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <cmath>
#include <vector>

ViewVector_d MIRCO::CreateMeshgrid(const int ngrid, const double GridSize)
{
  ViewVector_d meshgrid("meshgrid", ngrid);
  
  const double GridSize_2 = GridSize/2;
  Kokkos::parallel_for(
    "Create meshgrid", ngrid, KOKKOS_LAMBDA(const int i) {
      meshgrid(i) = GridSize_2 + i * GridSize;
    });
    
  return meshgrid;
}

MIRCO::TopologyMaxAndMean MIRCO::ComputeMaxAndMean(
    const Teuchos::SerialDenseMatrix<int, double>& topology)
{
  double zmax = -1.0;
  double zmean = 0.0;
#pragma omp parallel for schedule(guided, 16) reduction(+ : zmean) reduction(max : zmax)
  // Static and Guided seem even but Guided makes more sense
  for (int i = 0; i < topology.numRows(); i++)
  {
    for (int j = 0; j < topology.numCols(); j++)
    {
      zmean += topology(i, j);
      if (topology(i, j) > zmax)
      {
        zmax = topology(i, j);
      }
    }
  }

  zmean = zmean / (topology.numCols() * topology.numRows());
  return TopologyMaxAndMean{zmax, zmean};
}
