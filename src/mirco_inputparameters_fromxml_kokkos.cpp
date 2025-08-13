#include <pugixml.hpp>

#include "mirco_filesystem_utils.h"
#include "mirco_inputparameters_kokkos.h"

namespace
{
  double get_double(const pugi::xml_node& node, const std::string& tag)
  {
    auto ele = node.child(tag);
    if (!ele) throw std::runtime_error("Parameter \"" + tag + "\" not found");
    return std::stod(ele.attribute("value").value());
  }
  int get_int(const pugi::xml_node& node, const std::string& tag)
  {
    auto ele = node.child(tag);
    if (!ele) throw std::runtime_error("Parameter \"" + tag + "\" not found");
    return std::stoi(ele.attribute("value").value());
  }
  bool get_bool(const pugi::xml_node& node, const std::string& tag)
  {
    auto ele = node.child(tag);
    if (!ele) throw std::runtime_error("Parameter \"" + tag + "\" not found");
    std::string val = ele.attribute("value").value();
    return val == "true" || val == "1";
  }
  std::string get_string(const pugi::xml_node& node, const std::string& tag)
  {
    auto ele = node.child(tag);
    if (!ele) throw std::runtime_error("Parameter \"" + tag + "\" not found");
    return ele.attribute("value").value();
  }
}  // namespace

MIRCO::InputParameters::InputParameters(const std::string& inputFileName)
{
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file(inputFileName.c_str());
  if (!result)
    throw std::runtime_error(
        "Failed to parse XML input file: " + std::string(result.description()));

  pugi::xml_node root = doc.child("mirco_input");

  pugi::xml_node parameters = root.child("parameters");
  pugi::xml_node geoParams = parameters.child("geometrical_parameters");
  pugi::xml_node matParams = parameters.child("material_parameters");

  // Set the surface generator based on RandomTopologyFlag
  if (get_bool(root, "RandomTopologyFlag"))
  {
    *this = InputParameters(get_double(matParams, "E1"), get_double(matParams, "E2"),
        get_double(matParams, "nu1"), get_double(matParams, "nu2"),
        get_double(geoParams, "Tolerance"), get_double(geoParams, "Delta"),
        get_double(geoParams, "LateralLength"), get_int(geoParams, "Resolution"),
        get_double(geoParams, "InitialTopologyStdDeviation"),
        get_double(geoParams, "HurstExponent"), get_bool(root, "RandomSeedFlag"),
        get_int(root, "RandomGeneratorSeed"), get_int(root, "MaxIteration"),
        get_bool(root, "WarmStartingFlag"), get_bool(root, "PressureGreenFunFlag"));
  }
  else
  {
    std::string topology_file_path = get_string(root, "TopologyFilePath");
    // The following function generates the actual path of the topology file
    MIRCO::UTILS::ChangeRelativePath(topology_file_path, inputFileName);

    *this = InputParameters(get_double(matParams, "E1"), get_double(matParams, "E2"),
        get_double(matParams, "nu1"), get_double(matParams, "nu2"),
        get_double(geoParams, "Tolerance"), get_double(geoParams, "Delta"),
        get_double(geoParams, "LateralLength"), topology_file_path, get_int(root, "MaxIteration"),
        get_bool(root, "WarmStartingFlag"), get_bool(root, "PressureGreenFunFlag"));
  }
}
