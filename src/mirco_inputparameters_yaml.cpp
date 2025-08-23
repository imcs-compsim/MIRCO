#include <fstream>
#include <sstream>

#include "mirco_inputparameters_kokkos.h"
#include "mirco_utils_kokkos.h"

MIRCO::InputParameters::InputParameters(const std::string& inputFileName)
{
  std::ifstream fin(inputFileName);
  if (!fin) throw std::runtime_error("Cannot open file: " + inputFileName);

  std::stringstream ss;
  ss << fin.rdbuf();
  std::string inString = ss.str();

  ryml::Tree tree = ryml::parse_in_arena(c4::to_csubstr(inString));
  ryml::ConstNodeRef root = tree["mirco_input"];
  ryml::ConstNodeRef parameters = root["parameters"];
  ryml::ConstNodeRef geoParams = parameters["geometrical_parameters"];
  ryml::ConstNodeRef matParams = parameters["material_parameters"];

  if (root.invalid() || parameters.invalid() || geoParams.invalid() || matParams.invalid())
    throw std::runtime_error("Input incomplete");

  // Set the surface generator based on RandomTopologyFlag
  if (Utils::get_bool(root, "RandomTopologyFlag"))
  {
    *this = InputParameters(Utils::get_double(matParams, "E1"), Utils::get_double(matParams, "E2"),
        Utils::get_double(matParams, "nu1"), Utils::get_double(matParams, "nu2"),
        Utils::get_double(geoParams, "Tolerance"), Utils::get_double(geoParams, "Delta"),
        Utils::get_double(geoParams, "LateralLength"), Utils::get_int(geoParams, "Resolution"),
        Utils::get_double(geoParams, "InitialTopologyStdDeviation"),
        Utils::get_double(geoParams, "HurstExponent"), Utils::get_bool(root, "RandomSeedFlag"),
        Utils::get_int(root, "RandomGeneratorSeed"), Utils::get_int(root, "MaxIteration"),
        Utils::get_bool(root, "WarmStartingFlag"), Utils::get_bool(root, "PressureGreenFunFlag"));
  }
  else
  {
    std::string topology_file_path = Utils::get_string(root, "TopologyFilePath");
    // The following function generates the actual path of the topology file
    MIRCO::Utils::changeRelativePath(topology_file_path, inputFileName);

    *this = InputParameters(Utils::get_double(matParams, "E1"), Utils::get_double(matParams, "E2"),
        Utils::get_double(matParams, "nu1"), Utils::get_double(matParams, "nu2"),
        Utils::get_double(geoParams, "Tolerance"), Utils::get_double(geoParams, "Delta"),
        Utils::get_double(geoParams, "LateralLength"), topology_file_path,
        Utils::get_int(root, "MaxIteration"), Utils::get_bool(root, "WarmStartingFlag"),
        Utils::get_bool(root, "PressureGreenFunFlag"));
  }
}
