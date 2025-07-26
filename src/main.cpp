#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
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

  auto& topology = *(inputParams.topology_);
  auto max_and_mean = MIRCO::ComputeMaxAndMean(topology);

  // Initialise Pressure
  double pressure = 0.0;
  MIRCO::Evaluate(inputParams, pressure, max_and_mean.max_, meshgrid);

  std::cout << "Mean pressure is: " << std::to_string(pressure) << std::endl;

  const auto finish = std::chrono::high_resolution_clock::now();
  const double elapsedTime =
      std::chrono::duration_cast<std::chrono::duration<double>>(finish - start).count();
  std::cout << "Elapsed time is: " + std::to_string(elapsedTime) + "s." << std::endl;

  // Test for correct output if the result_description is given in the input file
  Teuchos::RCP<Teuchos::ParameterList> parameterList = Teuchos::rcp(new Teuchos::ParameterList());
  Teuchos::updateParametersFromXmlFile(inputFileName, parameterList.ptr());
  if (parameterList->isSublist("result_description"))
  {
    Teuchos::ParameterList& result_description = parameterList->sublist("result_description");
    double ExpectedPressure = result_description.get<double>("ExpectedPressure");
    double ExpectedPressureTolerance = result_description.get<double>("ExpectedPressureTolerance");
    if (std::abs(pressure - ExpectedPressure) > ExpectedPressureTolerance)
    {
      std::cerr << "The output pressure is incorrect" << std::endl;
      return EXIT_FAILURE;
    }
  }
}
