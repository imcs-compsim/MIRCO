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

  bool WarmStartingFlag;
  int Resolution;
  double MaxTopologyHeight;
  double nu1, nu2, CompositeYoungs, alpha, ElasticComplianceCorrection, GridSize, Tolerance, E1, E2,
      LateralLength, Delta;
  bool RandomTopologyFlag;
  bool RandomSeedFlag;
  double Hurst;
  std::string TopologyFilePath;
  int RandomGeneratorSeed;
  int MaxIteration;

  MIRCO::SetParameters(E1, E2, LateralLength, nu1, nu2, CompositeYoungs, alpha,
      ElasticComplianceCorrection, GridSize, Tolerance, Delta, TopologyFilePath, Resolution,
      MaxTopologyHeight, inputFileName, RandomTopologyFlag, Hurst, RandomSeedFlag,
      RandomGeneratorSeed, WarmStartingFlag, MaxIteration);

  // Identical Vectors/Matricies, therefore only created one here.
  int ngrid = int(ceil((LateralLength - (GridSize / 2)) / GridSize));
  std::vector<double> meshgrid(ngrid);
  MIRCO::CreateMeshgrid(meshgrid, ngrid, GridSize);

  // Setup Topology
  Epetra_SerialDenseMatrix topology;
  int N = pow(2, Resolution);
  topology.Shape(N + 1, N + 1);

  Teuchos::RCP<MIRCO::TopologyGeneration> surfacegenerator;
  // creating the correct surface object
  MIRCO::CreateSurfaceObject(Resolution, MaxTopologyHeight, Hurst, RandomSeedFlag, TopologyFilePath,
      RandomTopologyFlag, RandomGeneratorSeed, surfacegenerator);

  surfacegenerator->GetSurface(topology, MaxTopologyHeight);

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
