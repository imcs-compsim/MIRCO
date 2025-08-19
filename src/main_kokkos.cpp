#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "mirco_evaluate_kokkos.h"
#include "mirco_inputparameters_kokkos.h"
#include "mirco_kokkostypes_kokkos.h"
#include "mirco_topologyutilities_kokkos.h"
#include "mirco_utils_kokkos.h"

using namespace MIRCO;

int main(int argc, char* argv[])
{
  Kokkos::initialize(argc, argv);
  {
    std::cout << "-- Kokkos info --\n";
    std::cout << "threads_in_use = " << ExecSpace_Default_t::concurrency() << "\n";
    std::cout << "ExecSpace_Default = " << typeid(ExecSpace_Default_t).name() << "\n";
    std::cout << "ExecSpace_DefaultHost = " << typeid(ExecSpace_DefaultHost_t).name() << "\n";
    std::cout << "MemorySpace_Host_t = " << typeid(MemorySpace_Host_t).name() << "\n";
    std::cout << "MemorySpace_ofDefaultExec_t = " << typeid(MemorySpace_ofDefaultExec_t).name()
              << "\n\n";

    if (argc != 2) std::runtime_error("The code expects (only) an input file as argument");
    // Read the input file name from the command line
    std::string inputFileName = argv[1];

    const auto start = std::chrono::high_resolution_clock::now();

    InputParameters inputParams(inputFileName);

    ViewVector_d meshgrid_d = CreateMeshgrid(inputParams.N, inputParams.grid_size);
    auto maxAndMean = ComputeMaxAndMean(inputParams.topology_d);

    // Main evaluation agorithm
    double meanPressure;
    Evaluate(meanPressure, inputParams, maxAndMean.max, meshgrid_d);

    std::cout << "Mean pressure is: " << std::to_string(meanPressure) << std::endl;

    const auto finish = std::chrono::high_resolution_clock::now();
    const double elapsedTime =
        std::chrono::duration_cast<std::chrono::duration<double>>(finish - start).count();
    std::cout << "Elapsed time is: " + std::to_string(elapsedTime) + "s" << std::endl;

    // Test for correct output if the result_description is given in the input file
    {
      std::ifstream fin(inputFileName);
      if (!fin) throw std::runtime_error("Cannot open file: " + inputFileName);

      std::stringstream ss;
      ss << fin.rdbuf();
      std::string inString = ss.str();

      ryml::Tree tree = ryml::parse_in_arena(c4::to_csubstr(inString));
      ryml::ConstNodeRef root = tree["mirco_input"];
      ryml::ConstNodeRef resultDescription = root["result_description"];
      if (!resultDescription.invalid())
      {
        const double ExpectedPressure = Utils::get_double(resultDescription, "ExpectedPressure");
        const double ExpectedPressureTolerance =
            Utils::get_double(resultDescription, "ExpectedPressureTolerance");
        if (std::abs(meanPressure - ExpectedPressure) > ExpectedPressureTolerance)
        {
          std::cerr << "The output pressure does not match the expected result" << std::endl;
          return EXIT_FAILURE;
        }
        std::cout << "meanPressure=" << meanPressure << "\n";
        std::cout << "\tExpectedPressure=" << ExpectedPressure << "\n";
        std::cout << "\tExpectedPressureTolerance=" << ExpectedPressureTolerance << "\n";
      }
    }
  }
  Kokkos::finalize();
}
