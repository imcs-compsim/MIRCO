#include "mirco_exportvisualization.h"

#include "mirco_exportvisualization_vtk.h"

namespace MIRCO
{
  void ExportVisualization(const std::string& path, float gridSize, const ViewVectorInt_d activeSet,
      const std::vector<ViewMatrix_d>& otherFields, const std::vector<std::string>& otherFieldNames)
  {
    const int n = otherFields[0].extent(0);
    const int n2 = n * n;

    ViewVectorInt_h activeSet_h =
        Kokkos::create_mirror_view_and_copy(ExecSpace_DefaultHost_t(), activeSet);

    std::vector<unsigned char> activeSetPlain(n2, 0);

    for (int indA = 0; indA < activeSet_h.extent(0); ++indA) activeSetPlain[activeSet_h(indA)] = 1;

    std::vector<std::vector<float>> fieldsPlain(otherFields.size());

    for (int f = 0; f < otherFields.size(); ++f)
    {
      ViewMatrix_h field_h =
          Kokkos::create_mirror_view_and_copy(ExecSpace_DefaultHost_t(), otherFields[f]);

      fieldsPlain[f].resize(n2);

      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) fieldsPlain[f][i + n * j] = static_cast<float>(field_h(i, j));
    }

    ExportVisualizationVTK(path, gridSize, n, activeSetPlain, fieldsPlain, otherFieldNames);
  }
}  // namespace MIRCO
