#include "mirco_exportvisualization_vtk.h"

#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkZLibDataCompressor.h>

#include <iostream>

namespace MIRCO
{
  void ExportVisualizationVTK(const std::string& path, float gridSize, int n,
      const std::vector<unsigned char>& activeSet,
      const std::vector<std::vector<float>>& otherFields,
      const std::vector<std::string>& otherFieldNames)
  {
    const int n2 = n * n;

    vtkNew<vtkImageData> img;
    img->SetDimensions(n, n, 1);
    img->SetSpacing(gridSize, gridSize, 1.0);
    img->SetOrigin(0.0, 0.0, 0.0);

    // Active set
    {
      vtkNew<vtkUnsignedCharArray> vtkArr;
      vtkArr->SetName("Active Set");
      vtkArr->SetNumberOfComponents(1);
      vtkArr->SetNumberOfTuples(n2);

      for (int i = 0; i < n2; ++i) vtkArr->SetValue(i, activeSet[i]);

      img->GetPointData()->AddArray(vtkArr);
    }

    // Other fields
    for (int f = 0; f < otherFields.size(); ++f)
    {
      vtkNew<vtkFloatArray> vtkArr;
      vtkArr->SetName(otherFieldNames[f].c_str());
      vtkArr->SetNumberOfComponents(1);
      vtkArr->SetNumberOfTuples(n2);

      for (int i = 0; i < n2; ++i) vtkArr->SetValue(i, otherFields[f][i]);

      img->GetPointData()->AddArray(vtkArr);
    }

    vtkNew<vtkXMLImageDataWriter> w;
    w->SetFileName((path + ".vti").c_str());
    w->SetInputData(img);

    vtkNew<vtkZLibDataCompressor> z;
    w->SetCompressor(z);
    w->SetDataModeToAppended();
    w->EncodeAppendedDataOff();

    if (!w->Write()) std::cerr << "WARNING: ExportVisualization() failed to write to file\n\n";
  }

}  // namespace MIRCO
