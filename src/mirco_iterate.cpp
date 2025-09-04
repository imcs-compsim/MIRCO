#include "mirco_evaluate.h"
#include "mirco_iterate.h"


void MIRCO::Iterate(double& contactarea, double targetpressure, double initialguessDelta, double LateralLength, double GridSize,
    double Tolerance, int MaxIteration, double CompositeYoungs, bool WarmStartingFlag,
    double ElasticComplianceCorrection, Teuchos::SerialDenseMatrix<int, double>& topology,
    double zmax, std::vector<double>& meshgrid, bool PressureGreenFunFlag)
    {
        std::cout << "Starting the iteration to find the gap for a target pressure..." << std::endl;
        // //Pressure to solve for (hardcoded for testing)
        // double targetpressure = 0.000492321316110599;
        // //Initial guess for the gap (hardcoded for testing)
        // double InitialGuess = 15.90;

        double error = std::numeric_limits<double>::max();

        double iteration = 0;

        double CurrentDelta = initialguessDelta;
        double max_iter = 50000;

        while(error > 10e-6 && iteration < max_iter)
        {
            std::cout << "Current Delta: " << CurrentDelta << " micrometers." << std::endl;
            std::cout << "Current Iteration: " << iteration << std::endl;
            double OldDelta = CurrentDelta;

            double pressure = 0.0;
            MIRCO::Evaluate(contactarea, pressure, CurrentDelta, LateralLength, GridSize, Tolerance, MaxIteration,
               CompositeYoungs, WarmStartingFlag, ElasticComplianceCorrection, topology, zmax,
               meshgrid, PressureGreenFunFlag);
            std::cout << "Current Pressure: " << pressure << std::endl;

            double slightly_lower_pressure = 0;
            MIRCO::Evaluate(contactarea, slightly_lower_pressure, 0.999*CurrentDelta, LateralLength, GridSize, Tolerance, MaxIteration,
               CompositeYoungs, WarmStartingFlag, ElasticComplianceCorrection, topology, zmax,
               meshgrid, PressureGreenFunFlag);
            std::cout << "Slightly lower Pressure: " << slightly_lower_pressure << std::endl;

            CurrentDelta = CurrentDelta - ((pressure - targetpressure)*0.001*CurrentDelta)/(pressure - slightly_lower_pressure);
            error = std::abs((CurrentDelta - OldDelta)/OldDelta);
            std::cout << "Current Error: " << error << std::endl;
            std::cout << "----------------------------------------" << std::endl;
            iteration = iteration + 1;
        }

        std::cout << "Converged in " << iteration << " iterations." << std::endl;
        std::cout << "Contact area is: " << std::to_string(contactarea) << std::endl;

        

    }