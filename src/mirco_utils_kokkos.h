#ifndef SRC_UTILS_KOKKOS_H_
#define SRC_UTILS_KOKKOS_H_

#include <ryml.hpp>
#include <ryml_std.hpp>
#include <string>


namespace MIRCO
{
  namespace Utils
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
    void changeRelativePath(std::string& targetfilename, const std::string& sourcefilename);

    /*
     * \brief The following 4 functions are used to get the value of a parameter in a list (node) in
     * a yaml tree
     */
    std::string get_string(ryml::ConstNodeRef node, const std::string& key);
    bool get_bool(ryml::ConstNodeRef node, const std::string& key);
    double get_double(ryml::ConstNodeRef node, const std::string& key);
    int get_int(ryml::ConstNodeRef node, const std::string& key);
  }  // namespace Utils
}  // namespace MIRCO

#endif /* SRC_UTILS_KOKKOS_H_ */
