#ifndef SRC_FILESYSTEM_UTILS_H_
#define SRC_FILESYSTEM_UTILS_H_

#include <string>

namespace MIRCO
{
  namespace UTILS
  {
    /*
     * \brief Create a file path relative to the current directory by concatenating
     * two relative paths
     *
     * This function takes a sourcefilename to a file relative to the current
     * directory, and a targetfilename, that is relative to the sourcefilepath and
     * concatenates them. If the targetfilename is already an absolute path, it will
     * do nothing.
     */
    void ChangeRelativePath(std::string& targetfilename, const std::string& sourcefilename);
  }  // namespace UTILS
}  // namespace MIRCO

#endif /* SRC_FILESYSTEM_UTILS_H_ */
