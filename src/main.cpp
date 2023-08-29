#include <Teuchos_TestForException.hpp>
#include <chrono>
#include <iostream>
#include <string>

#include "mirco_evaluate.h"
#include "mirco_filesystem_utils.h"
#include "mirco_postprocess.h"
#include "mirco_setparameters.h"
#include "mirco_topology.h"
#include "mirco_topologyutilities.h"

int main(int argc, char* argv[])
{
  TEUCHOS_TEST_FOR_EXCEPTION(argc != 3, std::invalid_argument,
      "The code expects input file with a full path and a prefix for output files.");
  // reading the input file name from the command line
  std::string inputFileName = argv[1];
  std::string outputFileName = argv[2];

  // following function generates the actual path of the output file.
  UTILS::ChangeRelativePath(outputFileName, inputFileName);

  bool WarmStartingFlag = false;
  int Resolution = 0;
  double nu1 = 0.0, nu2 = 0.0, CompositeYoungs = 0.0, alpha = 0.0,
         ElasticComplianceCorrection = 0.0, GridSize = 0.0, Tolerance = 0.0, E1 = 0.0, E2 = 0.0,
         LateralLength = 0.0, Delta = 0.0;
  bool RandomTopologyFlag = false;
  bool RandomSeedFlag = false;
  double Hurst = 0.0;
  double InitialTopologyStdDeviation = 0.0;
  std::string TopologyFilePath = "";
  int RandomGeneratorSeed = 0;
  int MaxIteration = 0;

  MIRCO::SetParameters(E1, E2, LateralLength, nu1, nu2, CompositeYoungs, alpha,
      ElasticComplianceCorrection, GridSize, Tolerance, Delta, TopologyFilePath, Resolution,
      InitialTopologyStdDeviation, inputFileName, RandomTopologyFlag, Hurst, RandomSeedFlag,
      RandomGeneratorSeed, WarmStartingFlag, MaxIteration);

  // Identical Vectors/Matricies, therefore only created one here.
  int ngrid = int(ceil((LateralLength - (GridSize / 2)) / GridSize));
  std::vector<double> meshgrid(ngrid);
  MIRCO::CreateMeshgrid(meshgrid, ngrid, GridSize);

  // Setup Topology
  Teuchos::SerialDenseMatrix<int, double> topology;
  int N = pow(2, Resolution);
  topology.shape(N + 1, N + 1);

  Teuchos::RCP<MIRCO::TopologyGeneration> surfacegenerator;
  // creating the correct surface object
  MIRCO::CreateSurfaceObject(Resolution, InitialTopologyStdDeviation, Hurst, RandomSeedFlag,
      TopologyFilePath, RandomTopologyFlag, RandomGeneratorSeed, surfacegenerator);

  surfacegenerator->GetSurface(topology);

  auto max_and_mean = MIRCO::ComputeMaxAndMean(topology);

  // Number of contact points
  int nf = 0;
  // Coordinates of the points in contact in the last iteration.
  std::vector<double> xvf, yvf;
  // Contact force at (xvf,yvf) calculated in the last iteration.
  std::vector<double> pf;

  // Set start time
  const auto start = std::chrono::high_resolution_clock::now();

  MIRCO::Evaluate(Delta, LateralLength, GridSize, Tolerance, MaxIteration, CompositeYoungs,
      WarmStartingFlag, ElasticComplianceCorrection, topology, max_and_mean.max_, meshgrid, xvf,
      yvf, pf, nf);

  // Total MIRCO simulation time
  const auto finish = std::chrono::high_resolution_clock::now();
  const double elapsedTime =
      std::chrono::duration_cast<std::chrono::seconds>(finish - start).count();

  MIRCO::PostProcess(
      xvf, yvf, pf, nf, GridSize, ngrid, meshgrid, LateralLength, elapsedTime, outputFileName);
}
