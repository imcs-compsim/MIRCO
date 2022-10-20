#ifndef SRC_TOPOLOGYUTILITIES_H_
#define SRC_TOPOLOGYUTILITIES_H_

#include <Epetra_SerialSymDenseMatrix.h>
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
   * @param delta Grid size
   */
  void CreateMeshgrid(std::vector<double>& meshgrid, int ngrid, double delta);

  // forward declaration
  class TopologyGeneration;

  /**
   * @brief Create a Surface Object
   *
   * The aim of this this function is to create a surface object of correct object class depending
   * on the rmg_flag.
   *
   * @param resolution Resolution parameter
   * @param user_zmax Maximum height of the topology
   * @param Hurst Hurst Exponent (Used in random mid-point generator)
   * @param rand_seed_flag Set `true` to fix the seed to generate psuedo random topology to
   * reproduce results. Set `false` to use random seed.
   * @param zfilePath Path of the input file containing the topology relative to the build
   * directory.
   * @param rmg_flag Set `true` to use the Random-Midpoint Generator to generate the topology. Set
   * `false` to read topology from a file.
   * @param rmg_seed Seed for the random mid-point generator
   * @param surfacegenerator Surface object
   */
  void CreateSurfaceObject(int resolution, double& user_zmax, double Hurst, bool rand_seed_flag,
      std::string zfilePath, bool rmg_flag, int rmg_seed,
      std::shared_ptr<TopologyGeneration>& surfacegenerator);

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