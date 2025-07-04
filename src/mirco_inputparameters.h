#ifndef SRC_INPUTPARAMETERS_H_
#define SRC_INPUTPARAMETERS_H_

#include <Teuchos_ParameterList.hpp>  //# these are for the second ctor. though if we should move away from teuchos, we can change the SetParameters signature
#include <Teuchos_RCP.hpp>            //#
#include <Teuchos_XMLParameterListHelpers.hpp>  //#
#include <memory>
#include <string>

namespace MIRCO
{
  /**
   * @brief This class stores the input parameters
   *
   */
  class InputParameters
  {
   public:
    // ctor .xml (RMG or .dat)
    InputParameters(const std::string& inputFileName);
    // ctor no .xml; RMG
    InputParameters(double E1, double E2, double LateralLength, double nu1, double nu2,
        double& GridSize, double& Tolerance, double& Delta, int& Resolution,
        double& InitialTopologyStdDeviation, double& Hurst, bool& RandomSeedFlag,
        int& RandomGeneratorSeed, bool& WarmStartingFlag, int& MaxIteration,
        bool& PressureGreenFunFlag);
    // ctor no .xml; .dat
    InputParameters(double E1, double E2, double LateralLength, double nu1, double nu2,
        double& GridSize, double& Tolerance, double& Delta, std::string& TopologyFilePath,
        bool& WarmStartingFlag, int& MaxIteration, bool& PressureGreenFunFlag);



   private:
    void SetParameters(Teuchos::RCP<Teuchos::ParameterList> parameterList);

    double composite_youngs_ = 0.0, elastic_compliance_correction_ = 0.0, shape_factor_ = 0.0,
           grid_size_ = 0.0, lateral_length_ = 0.0, tolerance_ = 0.0, delta_ = 0.0;

    int max_iteration_ = 0;
    bool warm_starting_flag_ = false;
    bool PressureGreenFunFlag_ = false;
  };
}  // namespace MIRCO



// pass into MIRCO::Evaluate() just this object instead of all the params again

#endif  // SRC_INPUTPARAMETERS_H_
