#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
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
  double nu1 = 0.0, nu2 = 0.0, CompositeYoungs = 0.0, ShapeFactor = 0.0,
         ElasticComplianceCorrection = 0.0, GridSize = 0.0, Tolerance = 0.0, E1 = 0.0, E2 = 0.0,
         LateralLength = 0.0, Delta = 0.0;
  bool RandomTopologyFlag = false;
  bool RandomSeedFlag = false;
  double Hurst = 0.0;
  double InitialTopologyStdDeviation = 0.0;
  std::string TopologyFilePath = "";
  int RandomGeneratorSeed = 0;
  int MaxIteration = 0;
  bool PressureGreenFunFlag = false;

  double ExpectedPressure = -1.0;
  double ExpectedPressureTolerance = -1.0;

  Teuchos::RCP<Teuchos::ParameterList> parameterList = Teuchos::rcp(new Teuchos::ParameterList());
  Teuchos::updateParametersFromXmlFile(inputFileName, parameterList.ptr());

  MIRCO::SetParameters(E1, E2, LateralLength, nu1, nu2, CompositeYoungs, ShapeFactor,
      ElasticComplianceCorrection, GridSize, Tolerance, Delta, TopologyFilePath, Resolution,
      InitialTopologyStdDeviation, inputFileName, RandomTopologyFlag, Hurst, RandomSeedFlag,
      RandomGeneratorSeed, WarmStartingFlag, MaxIteration, PressureGreenFunFlag, ExpectedPressure,
      ExpectedPressureTolerance);

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

  // Initialise Pressure
  double pressure = 0.0;

  MIRCO::Evaluate(pressure, Delta, LateralLength, GridSize, Tolerance, MaxIteration,
      CompositeYoungs, WarmStartingFlag, ElasticComplianceCorrection, topology, max_and_mean.max_,
      meshgrid, PressureGreenFunFlag);

  std::cout << "Mean pressure is: " << std::to_string(pressure) << std::endl;

  const auto finish = std::chrono::high_resolution_clock::now();
  const double elapsedTime =
      std::chrono::duration_cast<std::chrono::duration<double>>(finish - start).count();
  std::cout << "Elapsed time is: " + std::to_string(elapsedTime) + "s." << std::endl;

  // Test for correct output if the result_description is given in the input file
  if (parameterList->isSublist("result_description") &&
      std::abs(pressure - ExpectedPressure) > ExpectedPressureTolerance)
  {
    std::cerr << "The output pressure is incorrect" << std::endl;
    return EXIT_FAILURE;
  }
}
