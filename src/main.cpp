#include <Teuchos_TestForException.hpp>
#include <chrono>
#include <iostream>
#include <string>

#include "mirco_evaluate.h"
#include "mirco_inputparameters.h"


int main(int argc, char* argv[])
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      argc != 2, std::invalid_argument, "The code expects (only) an input file as argument");
  // reading the input file name from the command line
  std::string inputFileName = argv[1];

  const auto start = std::chrono::high_resolution_clock::now();

  // Identical Vectors/Matricies, therefore only created one here.
  int ngrid =
      int(ceil((LateralLength - (GridSize / 2)) /
               GridSize));  // # wtf is this autism it is always something.5 and then we ceil
  std::cout << "ngrid=" << ngrid << "\n";
  MIRCO::InputParameters inputParams(inputFileName);
  auto topology = inputParams.topology_;
  std::vector<double> meshgrid(topology.numRows());
  MIRCO::CreateMeshgrid(meshgrid, ngrid, GridSize);  // #

  auto max_and_mean = MIRCO::ComputeMaxAndMean(topology);

  // Initialise Pressure
  double pressure = 0.0;

  MIRCO::Evaluate(pressure, Delta, LateralLength, GridSize, Tolerance, MaxIteration,
      CompositeYoungs, WarmStartingFlag, ElasticComplianceCorrection, topology, max_and_mean.max_,
      meshgrid, PressureGreenFunFlag);  // # make overloard constructor: either pass all these
                                        // params in or pass in an InputParameters object as well as
                                        // the necessary other params

  std::cout << "Mean pressure is: " << std::to_string(pressure) << std::endl;

  const auto finish = std::chrono::high_resolution_clock::now();
  const double elapsedTime =
      std::chrono::duration_cast<std::chrono::seconds>(finish - start).count();
  std::cout << "Elapsed time is: " + std::to_string(elapsedTime) + "s." << std::endl;
}
