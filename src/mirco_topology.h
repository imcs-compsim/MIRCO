#ifndef SRC_TOPOLOGY_H_
#define SRC_TOPOLOGY_H_

#include <Teuchos_SerialDenseMatrix.hpp>
#include <string>

namespace MIRCO
{
  /**
   * @brief Construct a Surface object
   *
   */
  class TopologyGeneration
  {
   public:
    int resolution;
    /**
     * @brief Get the Surface
     *
     * @param z Initialised topology matrix containing heights
     * @param zmax Maximum height of the topology
     */
    virtual void GetSurface(Teuchos::SerialDenseMatrix<int, double> &z) = 0;
    TopologyGeneration(int nn) { resolution = nn; }

    virtual ~TopologyGeneration() = default;
  };

  class ReadFile : public TopologyGeneration
  {
   public:
    std::string TopologyFilePath;

    void GetSurface(Teuchos::SerialDenseMatrix<int, double> &z) override;

    /**
     * @brief Construct a Surface object by reading topology from an input file.
     *
     * @param nn Resolution parameter
     * @param ffilepath Path of the input file containing the topology relative to the build
     * directory.
     */
    ReadFile(int nn, std::string ffilepath) : TopologyGeneration(nn)
    {
      TopologyFilePath = ffilepath;
    }
  };

  class Rmg : public TopologyGeneration
  {
   public:
    double InitialTopologyStdDeviation;
    double Hurst;
    bool RandomSeedFlag;
    int RandomGeneratorSeed;

    void GetSurface(Teuchos::SerialDenseMatrix<int, double> &z) override;

    /**
     * @brief Construct a Surface object using Random Midpoint Generator
     *
     * @param nn Resolution parameter
     * @param InStdDev Initial Standard deviation for the random-midpoint generator
     * [micrometers]
     * @param HH Hurst exponent
     * @param rsf Random Seed Flag
     * @param rmgs eed for the random mid-point generator
     */
    Rmg(int nn, double InStdDev, double HH, bool rsf, int rmgs) : TopologyGeneration(nn)
    {
      InitialTopologyStdDeviation = InStdDev;
      Hurst = HH;
      RandomSeedFlag = rsf;
      RandomGeneratorSeed = rmgs;
    }
  };
}  // namespace MIRCO

#endif  // SRC_TOPOLOGY_H_
