#include <Teuchos_TestForException.hpp>
#include <string>
#include "evaluate.h"
#include "setparameters.h"
#include "topology.h"
#include "topologyutilities.h"

int main(int argc, char *argv[])
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      argc != 2, std::invalid_argument, "The code expects (only) an input file as argument");
  // reading the input file name from the command line
  std::string inputFileName = argv[1];

  bool flagwarm;
  int resolution;
  double nu1, nu2, G1, G2, E, alpha, k_el, delta, nnodi, to1, E1, E2, lato, errf, Delta;
  bool rmg_flag;
  bool rand_seed_flag;
  double Hurst;
  std::string zfilePath;
  int rmg_seed;
  int max_iter;

  MIRCO::SetParameters(E1, E2, lato, nu1, nu2, G1, G2, E, alpha, k_el, delta, nnodi, errf, to1,
      Delta, zfilePath, resolution, inputFileName, rmg_flag, Hurst, rand_seed_flag, rmg_seed,
      flagwarm, max_iter);

  // Identical Vectors/Matricies, therefore only created one here.
  int iter = int(ceil((lato - (delta / 2)) / delta));
  std::vector<double> meshgrid(iter);
  MIRCO::CreateMeshgrid(meshgrid, iter, delta);

  // Setup Topology
  Epetra_SerialDenseMatrix topology;
  int N = pow(2, resolution);
  topology.Shape(N + 1, N + 1);

  std::shared_ptr<MIRCO::TopologyGeneration> surfacegenerator;
  // creating the correct surface object
  MIRCO::CreateSurfaceObject(
      resolution, Hurst, rand_seed_flag, zfilePath, rmg_flag, rmg_seed, surfacegenerator);

  surfacegenerator->GetSurface(topology);

  double zmax = 0;
  double zmean = 0;

  MIRCO::ComputeMaxAndMean(topology, zmax, zmean);

  double pressure = 0.0;

  MIRCO::Evaluate(pressure, Delta, lato, delta, errf, to1, max_iter, E, flagwarm, k_el, topology,
      zmax, meshgrid);
}
