#ifndef SRC_TOPOLOGYUTILITIES_H_
#define SRC_TOPOLOGYUTILITIES_H_

#include <Epetra_SerialSymDenseMatrix.h>
#include <Teuchos_RCP.hpp>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace MIRCO
{
  /**
   * @brief Create a Meshgrid vector of the surface
   *
   * Creates a meshgrid with coordinates. Since the vector is identical in both directions,
   * therefore only one is created.
   *
   * @param meshgrid Meshgrid vector
   * @param ngrid Number of grid points in one direction
   * @param GridSize Grid size
   */
  void CreateMeshgrid(std::vector<double>& meshgrid, int ngrid, double GridSize);

  // forward declaration
  class TopologyGeneration;

  /**
   * @brief Create a Surface Object
   *
   * The aim of this this function is to create a surface object of correct object class depending
   * on the RandomTopologyFlag.
   *
   * @param Resolution Resolution parameter
   * @param InitialStdDeviation Initial Standard deviation for the random-midpoint generator
   * [micrometers]
   * @param Hurst Hurst Exponent (Used in random mid-point generator)
   * @param RandomSeedFlag Set `true` to fix the seed to generate psuedo random topology to
   * reproduce results. Set `false` to use random seed.
   * @param TopologyFilePath Path of the input file containing the topology relative to the build
   * directory.
   * @param RandomTopologyFlag Set `true` to use the Random-Midpoint Generator to generate the
   * topology. Set `false` to read topology from a file.
   * @param RandomGeneratorSeed Seed for the random mid-point generator
   * @param surfacegenerator Surface object
   */
  void CreateSurfaceObject(int Resolution, double InitialStdDeviation, double Hurst,
      bool RandomSeedFlag, std::string TopologyFilePath, bool RandomTopologyFlag,
      int RandomGeneratorSeed, Teuchos::RCP<TopologyGeneration>& surfacegenerator);

  /**
   * @brief Compute the maximum height (zmax) and the mean height (zmean) of the topology.
   *
   * @param topology Topology matrix containing heights
   * @param zmax Maximum height
   * @param zmean Mean height
   */
  void ComputeMaxAndMean(Epetra_SerialDenseMatrix topology, double& zmax, double& zmean);
}  // namespace MIRCO

#endif  // SRC_TOPOLOGYUTILITIES_H_