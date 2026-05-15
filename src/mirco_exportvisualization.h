#ifndef SRC_EXPORTVISUALIZATION_H_
#define SRC_EXPORTVISUALIZATION_H_

#include <string>
#include <vector>

#include "mirco_kokkostypes.h"

namespace MIRCO
{
  void ExportVisualization(const std::string& path, float gridSize, const ViewVectorInt_d activeSet,
      const std::vector<ViewMatrix_d>& otherFields,
      const std::vector<std::string>& otherFieldNames);
}  // namespace MIRCO

#endif  // SRC_EXPORTVISUALIZATION_H_
