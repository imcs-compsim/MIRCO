#include "mirco_filesystem_utils.h"
#include <filesystem>

void UTILS::ChangeRelativePath(std::string& targetfilename, const std::string& sourcefilename)
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
