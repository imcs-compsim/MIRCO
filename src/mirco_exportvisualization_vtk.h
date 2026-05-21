#ifndef SRC_EXPORTVISUALIZATIONVTK_H_
#define SRC_EXPORTVISUALIZATIONVTK_H_

#include <string>
#include <vector>

namespace MIRCO
{
  void ExportVisualizationVTK(const std::string& path, float gridSize, int n,
      const std::vector<unsigned char>& activeSet,
      const std::vector<std::vector<float>>& otherFields,
      const std::vector<std::string>& otherFieldNames);
}  // namespace MIRCO

#endif  // SRC_EXPORTVISUALIZATIONVTK_H_
