#ifndef SRC_INPUTPARAMETERS_H_
#define SRC_INPUTPARAMETERS_H_

#include <Teuchos_ParameterList.hpp>  //# these are for the second ctor. though if we should move away from teuchos, we can change the SetParameters signature
#include <Teuchos_RCP.hpp>            //#
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>  //#
#include <memory>
#include <string>

namespace MIRCO
{
  /**
   * @brief This class stores the input parameters
   *
   */
  class InputParameters  // # at this point we should consider renaming to getTopology or something;
                         // but yeah, the whole goal of the input procedure is to get the topology
                         // anyways, so I think it is ok to have this all in a single class for
                         // conciseness and clarity--you won't ever be using anything else anyways;
                         // it is like impossible to separate it more while not having redundancies
  {
   public:
    // ctor .xml (RMG or .dat)
    InputParameters(const std::string& inputFileName);
    // ctor no .xml; RMG
    InputParameters(double E1, double E2, double nu1, double nu2, double Tolerance, double Delta,
        int Resolution, double LateralLength, double InitialTopologyStdDeviation, double Hurst,
        bool RandomSeedFlag, int RandomGeneratorSeed, int MaxIteration, bool WarmStartingFlag,
        bool PressureGreenFunFlag);
    // ctor no .xml; .dat
    InputParameters(double E1, double E2, double nu1, double nu2, double Tolerance, double Delta,
        double LateralLength, const std::string& TopologyFilePath, int MaxIteration,
        bool WarmStartingFlag, bool PressureGreenFunFlag);



    // In my opinion, we do not need getters
   public:
    double composite_youngs_ = 0.0, elastic_compliance_correction_ = 0.0, shape_factor_ = 0.0,
           tolerance_ = 0.0, delta_ = 0.0, lateral_length_ = 0.0, grid_size_ = 0.0;
    int max_iteration_ = 0;
    bool warm_starting_flag_ = false;
    bool pressure_greenfun_flag_ = false;

    Teuchos::SerialDenseMatrix<int, double> topology_;
  };
}  // namespace MIRCO



// pass into MIRCO::Evaluate() just this object instead of all the params again

#endif  // SRC_INPUTPARAMETERS_H_
