#ifndef SRC_TOPOLOGY_H_
#define SRC_TOPOLOGY_H_

#include <Epetra_SerialDenseMatrix.h>
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
    virtual void GetSurface(Epetra_SerialDenseMatrix &z, double &zmax) = 0;
    TopologyGeneration(int nn) { resolution = nn; }
  };

  class ReadFile : public TopologyGeneration
  {
   public:
    std::string TopologyFilePath;
    void GetSurface(Epetra_SerialDenseMatrix &z, double &zmax) override;
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
    double user_zmax;
    double Hurst;
    bool RandomSeedFlag;
    int RandomGeneratorSeed;
    void GetSurface(Epetra_SerialDenseMatrix &z, double &zmax) override;
    /**
     * @brief Construct a Surface object using Random Midpoint Generator
     *
     * @param nn Resolution parameter
     * @param u_zmax Maximum height of the topology
     * @param HH Hurst exponent
     * @param rsf Random Seed Flag
     * @param rmgs eed for the random mid-point generator
     */
    Rmg(int nn, double u_zmax, double HH, bool rsf, int rmgs) : TopologyGeneration(nn)
    {
      user_zmax = u_zmax;
      Hurst = HH;
      RandomSeedFlag = rsf;
      RandomGeneratorSeed = rmgs;
    }
  };
}  // namespace MIRCO

#endif  // SRC_TOPOLOGY_H_
