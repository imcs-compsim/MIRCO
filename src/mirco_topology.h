#ifndef SRC_TOPOLOGY_H_
#define SRC_TOPOLOGY_H_

#include <string>

#include "mirco_kokkostypes.h"

namespace MIRCO
{
  /**
   * @brief Construct a Surface object by reading topology from an input (.dat) file.
   *
   * @param[in] filepath Path of the input file containing the topology, relative to the calling
   * directory or absolute
   * @param[out] N Dimension of the surface (number of elements/cells/height nodes per sidelength)
   *
   * @return Topology heightfield matrix
   */
  ViewMatrix_h CreateSurfaceFromFile(const std::string& filepath, int& N);

  /**
   * @brief Construct a Surface object using Random Midpoint Generator
   *
   * @param[in] resolution Resolution parameter
   * @param[in] initialTopologyStdDeviation Initial Standard deviation for the random-midpoint
   * generator [micrometers]
   * @param[in] hurst Hurst exponent
   * @param[in] randomSeedFlag Random Seed Flag
   * @param[in] randomGeneratorSeed Seed for the random mid-point generator
   *
   * @return Topology heightfield matrix
   */
  ViewMatrix_h CreateRmgSurface(int resolution, double initialTopologyStdDeviation, double hurst,
      bool randomSeedFlag, int randomGeneratorSeed);

}  // namespace MIRCO

#endif  // SRC_TOPOLOGY_H_
