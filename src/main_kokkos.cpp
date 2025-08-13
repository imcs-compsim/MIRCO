#include <chrono>
#include <iostream>
#include <pugixml.hpp>
#include <string>

#include "mirco_evaluate_kokkos.h"
#include "mirco_inputparameters_kokkos.h"
#include "mirco_kokkostypes_kokkos.h"
#include "mirco_topologyutilities_kokkos.h"

int main(int argc, char* argv[])
{
  using namespace MIRCO;

  Kokkos::initialize(argc, argv);
  {
    int threads_in_use = ExecSpace_Default_t::concurrency();
    std::cout << "\nthreads_in_use=" << threads_in_use << "\n";
    std::cout << "\nExecSpace_Default= " << typeid(ExecSpace_Default_t).name() << "\n";
    std::cout << "\nExecSpace_DefaultHost= " << typeid(ExecSpace_DefaultHost_t).name() << "\n\n";
    std::cout << "\nMemorySpace_Host_t= " << typeid(MemorySpace_Host_t).name() << "\n";
    std::cout << "\nMemorySpace_ofDefaultExec_t= " << typeid(MemorySpace_ofDefaultExec_t).name()
              << "\n\n";

    if (argc != 2) std::runtime_error("The code expects (only) an input file as argument");
    // reading the input file name from the command line
    std::string inputFileName = argv[1];

    const auto start = std::chrono::high_resolution_clock::now();

    MIRCO::InputParameters inputParams(inputFileName);

    ViewVector_d meshgrid_d = MIRCO::CreateMeshgrid(inputParams.N, inputParams.grid_size);
    auto maxAndMean = MIRCO::ComputeMaxAndMean(inputParams.topology_d);

    // Main evaluation agorithm
    double meanPressure;
    MIRCO::Evaluate(meanPressure, inputParams, maxAndMean.max, meshgrid_d);

    std::cout << "Mean pressure is: " << std::to_string(meanPressure) << std::endl;

    const auto finish = std::chrono::high_resolution_clock::now();
    const double elapsedTime =
        std::chrono::duration_cast<std::chrono::duration<double>>(finish - start).count();
    std::cout << "Elapsed time is: " + std::to_string(elapsedTime) + "s." << std::endl;

    // Test for correct output if the result_description is given in the input file
    {
      pugi::xml_document doc;
      doc.load_file(inputFileName.c_str());
      auto result_description = doc.child("mirco_input").child("result_description");
      if (result_description)
      {
        const double ExpectedPressure =
            std::stod(result_description.child("ExpectedPressure").attribute("value").value());
        const double ExpectedPressureTolerance = std::stod(
            result_description.child("ExpectedPressureTolerance").attribute("value").value());
        if (std::abs(meanPressure - ExpectedPressure) > ExpectedPressureTolerance)
        {
          std::cerr << "The output pressure does not match the expected result" << std::endl;
          return EXIT_FAILURE;
        }
      }
    }
  }
  Kokkos::finalize();
}
