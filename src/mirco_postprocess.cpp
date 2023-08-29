#include "mirco_postprocess.h"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

long findClosestIndex(const std::vector<double>& sorted_array, double x)
{
  auto iter_geq = std::lower_bound(sorted_array.begin(), sorted_array.end(), x);

  if (iter_geq == sorted_array.begin())
  {
    return 0;
  }

  double a = *(iter_geq - 1);
  double b = *(iter_geq);

  if (fabs(x - a) < fabs(x - b))
  {
    return iter_geq - sorted_array.begin() - 1;
  }
  return iter_geq - sorted_array.begin();
}

void MIRCO::PostProcess(std::vector<double> xvf, std::vector<double> yvf, std::vector<double> pf,
    int nf, double GridSize, int ngrid, std::vector<double>& meshgrid, double LateralLength,
    double elapsedTime, std::string outputFileName)
{
  double forceArray[ngrid][ngrid] = {0};
  int ix;
  int iy;
  double areaContact;

  for (int i = 0; i < nf; i++)
  {
    ix = findClosestIndex(meshgrid, xvf[i]);
    iy = findClosestIndex(meshgrid, yvf[i]);
    forceArray[ix][ngrid-1-iy] = pf[i];
    // std::cout << ix << "\t" << iy  << std::endl;
    // std::cout << xvf[i] << "\t" << yvf[i] << "\t" << pf[i] << std::endl;
  }

  // Store force at every point (including contact and non-contact points)
  std::ofstream outputForceFile(outputFileName + "_force" + ".dat");

  for (int i = 0; i < ngrid; i++)
  {
    for (int j = 0; j < ngrid; j++)
    {
      outputForceFile << forceArray[j][ngrid - 1 - i] << ' ';
    }
    outputForceFile << "\n";
  }
  outputForceFile.close();

  // Calculate average pressure via total force
  auto totalForce = std::accumulate(pf.begin(), pf.end(), 0.0);
  double pressure = totalForce / pow(LateralLength, 2);

  // Total contact area
  areaContact = nf * (pow(GridSize, 2) / pow(LateralLength, 2)) * 100;

  // Store overall information
  std::ofstream infoFile(outputFileName + "_info" + ".csv");
  infoFile << "total_force"
           << "\t"
           << "total_area"
           << "\t"
           << "ave_pressure"
           << "\t"
           << "n_contact"
           << "\t"
           << "time"
           << "\n";
  infoFile << totalForce << "\t" << areaContact << "\t" << pressure << "\t" << nf << "\t"
           << std::to_string(elapsedTime);
  infoFile.close();
}