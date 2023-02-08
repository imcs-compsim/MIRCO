#include <Teuchos_TestForException.hpp>
#include <chrono>
#include <iostream>
#include <string>

#include "mirco_evaluate.h"
#include "mirco_setparameters.h"
#include "mirco_topology.h"
#include "mirco_topologyutilities.h"


int main(int argc, char* argv[])
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      argc != 2, std::invalid_argument, "The code expects (only) an input file as argument");
  // reading the input file name from the command line
  std::string inputFileName = argv[1];

  const auto start = std::chrono::high_resolution_clock::now();

  bool WarmStartingFlag = false;
  int Resolution = 0;
  double MaxTopologyHeight = 0.0;
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

  double MeanTopologyHeight = 0.0;
  MIRCO::ComputeMaxAndMean(topology, MaxTopologyHeight, MeanTopologyHeight);

  // Initialise Pressure
  double pressure = 0.0;

  MIRCO::Evaluate(pressure, Delta, LateralLength, GridSize, Tolerance, MaxIteration,
      CompositeYoungs, WarmStartingFlag, ElasticComplianceCorrection, topology, MaxTopologyHeight,
      meshgrid);

  std::cout << "Mean pressure is: " << std::to_string(pressure) << std::endl;

  const auto finish = std::chrono::high_resolution_clock::now();
  const double elapsedTime =
      std::chrono::duration_cast<std::chrono::seconds>(finish - start).count();
  std::cout << "Elapsed time is: " + std::to_string(elapsedTime) + "s." << std::endl;
}
