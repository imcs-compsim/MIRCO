#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "mirco_evaluate.h"
#include "mirco_inputparameters.h"
#include "mirco_kokkostypes.h"
#include "mirco_topologyutilities.h"
#include "mirco_utils.h"

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

    ViewVector_d meshgrid = CreateMeshgrid(inputParams.N, inputParams.grid_size);
    const auto maxAndMean = ComputeMaxAndMean(inputParams.topology);

    // Main evaluation agorithm
    double meanPressure, effectiveContactAreaFraction;
    Evaluate(meanPressure, effectiveContactAreaFraction, inputParams, maxAndMean.max, meshgrid);

    std::cout << std::setprecision(16) << "Mean pressure is: " << meanPressure
              << "\nEffective contact area fraction is: " << effectiveContactAreaFraction
              << std::endl;

    const auto finish = std::chrono::high_resolution_clock::now();
    const double elapsedTime =
        std::chrono::duration_cast<std::chrono::duration<double>>(finish - start).count();
    std::cout << "Elapsed time is: " + std::to_string(elapsedTime) + "s" << std::endl;

    // Test for correct output if the result_description is given in the input file
    {
      std::ifstream fin(inputFileName);
      if (!fin) throw std::runtime_error("Cannot open input file: " + inputFileName);

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
        const double ExpectedEffectiveContactAreaFraction =
            Utils::get_double(resultDescription, "ExpectedEffectiveContactAreaFraction");
        const double ExpectedEffectiveContactAreaFractionTolerance =
            Utils::get_double(resultDescription, "ExpectedEffectiveContactAreaFractionTolerance");

        if (std::abs(meanPressure - ExpectedPressure) > ExpectedPressureTolerance)
        {
          std::cerr << "The output pressure does not match the expected result." << std::endl;
          return EXIT_FAILURE;
        }
        if (std::abs(effectiveContactAreaFraction - ExpectedEffectiveContactAreaFraction) >
            ExpectedEffectiveContactAreaFractionTolerance)
        {
          std::cerr << "The output effective contact area does not match the expected result."
                    << std::endl;
          return EXIT_FAILURE;
        }
        std::cout << "meanPressure = " << meanPressure << "\n";
        std::cout << "\tExpectedPressure = " << ExpectedPressure << "\n";
        std::cout << "\tExpectedPressureTolerance = " << ExpectedPressureTolerance << "\n";
        std::cout << "effectiveContactArea = " << effectiveContactAreaFraction << "\n";
        std::cout << "\tExpectedEffectiveContactAreaFraction = "
                  << ExpectedEffectiveContactAreaFraction << "\n";
        std::cout << "\tExpectedEffectiveContactAreaFractionTolerance = "
                  << ExpectedEffectiveContactAreaFractionTolerance << std::endl;
      }
    }
  }
  Kokkos::finalize();
}
