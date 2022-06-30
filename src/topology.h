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
     */
    virtual void GetSurface(Epetra_SerialDenseMatrix &z) = 0;
    TopologyGeneration(int nn) { resolution = nn; }
  };

  class ReadFile : public TopologyGeneration
  {
   public:
    std::string filepath;
    void GetSurface(Epetra_SerialDenseMatrix &z) override;
    /**
     * @brief Construct a Surface object by reading topology from an input file.
     *
     * @param nn Resolution parameter
     * @param ffilepath Path of the input file containing the topology relative to the build
     * directory.
     */
    ReadFile(int nn, std::string ffilepath) : TopologyGeneration(nn) { filepath = ffilepath; }
  };

  class Rmg : public TopologyGeneration
  {
   public:
    double Hurst;
    bool rand_seed_flag;
    int rmg_seed;
    void GetSurface(Epetra_SerialDenseMatrix &z) override;
    /**
     * @brief Construct a Surface object using Random Midpoint Generator
     *
     * @param nn Resolution parameter
     * @param HH Hurst exponent
     * @param rsf Random Seed Flag
     * @param rmgs eed for the random mid-point generator
     */
    Rmg(int nn, double HH, bool rsf, int rmgs) : TopologyGeneration(nn)
    {
      Hurst = HH;
      rand_seed_flag = rsf;
      rmg_seed = rmgs;
    }
  };
}  // namespace MIRCO

#endif  // SRC_TOPOLOGY_H_
