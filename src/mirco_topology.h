#ifndef SRC_TOPOLOGY_H_
#define SRC_TOPOLOGY_H_

#include <Teuchos_SerialDenseMatrix.hpp>
#include <string>

namespace MIRCO
{

  /**
   * @brief Construct a Surface object by reading topology from an input (.dat) file.
   *
   * @param filepath Path of the input file containing the topology, relative to the calling
   * directory or absolute
   * @param[out] N Dimension of the surface (number of elements/cells/height nodes per sidelength)
   */
  Teuchos::SerialDenseMatrix<int, double> CreateSurfaceFromFile(
      const std::string& filepath, int& N);

  /**
   * @brief Construct a Surface object using Random Midpoint Generator
   *
   * @param resolution Resolution parameter
   * @param initialTopologyStdDeviation Initial Standard deviation for the random-midpoint generator
   * [micrometers]
   * @param hurst Hurst exponent
   * @param randomSeedFlag Random Seed Flag
   * @param randomGeneratorSeed Seed for the random mid-point generator
   */
  Teuchos::SerialDenseMatrix<int, double> CreateRmgSurface(int resolution,
      double initialTopologyStdDeviation, double hurst, bool randomSeedFlag,
      int randomGeneratorSeed);

}  // namespace MIRCO

#endif  // SRC_TOPOLOGY_H_
