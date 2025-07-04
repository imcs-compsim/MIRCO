#include <Teuchos_TestForException.hpp>
#include <chrono>
#include <iostream>
#include <string>

#include "mirco_evaluate.h"
#include "mirco_inputparameters.h"
#include "mirco_topologyutilities.h"


int main(int argc, char* argv[])
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      argc != 2, std::invalid_argument, "The code expects (only) an input file as argument");
  // reading the input file name from the command line
  std::string inputFileName = argv[1];

  const auto start = std::chrono::high_resolution_clock::now();

  MIRCO::InputParameters inputParams(inputFileName);

  // Identical Vectors/Matricies, therefore only created one here.
  auto meshgrid = MIRCO::CreateMeshgrid(inputParams.N_, inputParams.grid_size_);

  auto& topology = inputParams.topology_;
  auto max_and_mean = MIRCO::ComputeMaxAndMean(topology);

  // Initialise Pressure
  double pressure = 0.0;
  MIRCO::Evaluate(inputParams, pressure, max_and_mean.max_, meshgrid);

  std::cout << "Mean pressure is: " << std::to_string(pressure) << std::endl;

  const auto finish = std::chrono::high_resolution_clock::now();
  const double elapsedTime =
      std::chrono::duration_cast<std::chrono::seconds>(finish - start).count();
  std::cout << "Elapsed time is: " + std::to_string(elapsedTime) + "s." << std::endl;
}
