#include <filesystem>

#include "mirco_utils_kokkos.h"

void MIRCO::Utils::changeRelativePath(
    std::string& targetfilename, const std::string& sourcefilename)
{
  std::filesystem::path targetfilepath = targetfilename;

  if (targetfilepath.is_relative())
  {
    std::filesystem::path sourcefilepath = sourcefilename;
    std::filesystem::path root_dir = sourcefilepath.parent_path();
    root_dir += "/";
    root_dir += targetfilepath;
    targetfilename = root_dir.c_str();
  }
}

std::string MIRCO::Utils::get_string(ryml::ConstNodeRef node, const std::string& key)
{
  auto child = node[c4::to_csubstr(key)];
  if (child.invalid()) throw std::runtime_error("Parameter \"" + key + "\" not found");
  c4::csubstr v = child.val();
  return std::string(v.str, v.len);
}
bool MIRCO::Utils::get_bool(ryml::ConstNodeRef node, const std::string& key)
{
  std::string s = get_string(node, key);
  if (s == "true" || s == "True" || s == "1") return true;
  if (s == "false" || s == "False" || s == "0") return false;

  throw std::runtime_error("Parameter \"" + key + "\" has invalid bool value");
}
double MIRCO::Utils::get_double(ryml::ConstNodeRef node, const std::string& key)
{
  std::string s = get_string(node, key);
  return std::stod(s);
}
int MIRCO::Utils::get_int(ryml::ConstNodeRef node, const std::string& key)
{
  std::string s = get_string(node, key);
  return std::stoi(s);
}
