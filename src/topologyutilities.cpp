#include "topologyutilities.h"
#include <Epetra_SerialSymDenseMatrix.h>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include "topology.h"

void CreateMeshgrid(std::vector<double>& x, int iter, double delta)
{
  // the aim of this function is to create a mesh grid of the surface.
#pragma omp parallel for schedule(static, 16)  // Same amount of work -> static
  for (int i = 0; i < iter; i++)
  {
    x[i] = (delta / 2) + i * delta;
  }
}

void CreateSurfaceObject(int n, double Hurst, bool rand_seed_flag, std::string zfilePath,
    bool rmg_flag, int rmg_seed, std::shared_ptr<TopologyGeneration>& surfacegenerator)
{
  // The aim of this this function is to create a surface object of correct object class depending
  // on the rmg_flag.
  if (rmg_flag)
  {
    surfacegenerator = std::shared_ptr<Rmg>(new Rmg(n, Hurst, rand_seed_flag, rmg_seed));
  }
  else
  {
    surfacegenerator = std::shared_ptr<ReadFile>(new ReadFile(n, zfilePath));
  }
}

void ComputeMaxAndMean(Epetra_SerialDenseMatrix topology, double& zmax, double& zmean)
{
  // The aim of this function is to compute the maximum height (zmax) and the mean height (zmean) of
  // the topology.
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