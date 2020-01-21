#include "QuaC_Pulse_Visitor.hpp"
#include "xacc.hpp"
#include <chrono>
#include <iomanip>
extern "C" {
#include "interface_xacc_ir.h"
}

namespace {
   void writeTimesteppingDataToCsv(const std::string& in_fileName, const TSData* const in_tsData, int in_nbSteps)
   {
      if (in_nbSteps < 1)
      {
         return;
      }

      const auto stringEndsWith = [](const std::string& in_string, const std::string& in_ending) {
         if (in_ending.size() > in_string.size()) 
         {
            return false;
         }

         return std::equal(in_ending.rbegin(), in_ending.rend(), in_string.rbegin());   
      };

      const auto getCurrentTimeString = [](){
         const auto currentTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
         std::stringstream ss;
         ss << std::put_time(std::localtime(&currentTime), "%Y%m%d_%X");
         return ss.str();
      };

      std::ofstream outputFile;
      std::string fileName = stringEndsWith(in_fileName, ".csv") ?  in_fileName.substr(0, in_fileName.size() - 4) : in_fileName;
      // Add a timestamp to prevent duplicate
      fileName += "_";
      fileName += getCurrentTimeString();
      fileName += ".csv";
      outputFile.open(fileName, std::ofstream::out);

      if(!outputFile.is_open())
      {
         std::cout << "Cannot open CSV file '" << fileName << "' for writing!\n";
		   return;
	   }

      const auto nbPopulations = in_tsData[0].nbPops;
      // Header
      outputFile << "Time, ";      
      for(int j = 0; j < nbPopulations; ++j)
      {
         outputFile << "Population[" << j << "], ";
      }
      outputFile << "\n";

      // Data
      for (int i = 0; i < in_nbSteps; ++i)
      {
         const TSData dataAtStep = in_tsData[i];
         outputFile << dataAtStep.time << ", ";
         for(int j = 0; j < nbPopulations; ++j)
         {
            outputFile << dataAtStep.populations[j] << ", ";
         }
         outputFile << "\n";
      }
      outputFile.close();
      std::cout << "Time-stepping data is written to file '" << fileName << "'\n";
   }        
}

// Export the time-stepping data to csv (e.g. for plotting)
// Uncomment to get the data exported
// #define EXPORT_TS_DATA_AS_CSV

namespace QuaC {
   void PulseVisitor::initialize(std::shared_ptr<AcceleratorBuffer> buffer, const HeterogeneousMap& in_params) 
   {
      // Debug
      std::cout << "Initialize Pulse simulator \n";
         
      // Initialize some params for testing
      TODO(Get rid of these default params and enforce that they must be set upstream.)
      double dt = 0.01;
      double stopTime = 8.0;
      int stepMax = 100000000;
      double nu = 5.0;
      double omega = 2 * M_PI * nu;
      // Qubit decay: just use a very small value
      double kappa = 0.0001;

      // Create a pulse controller
      m_pulseChannelController = std::make_unique<PulseChannelController>();
      {
         // Initialize the pulse constroller:
         // Note: the data will eventually be loaded from XACC. For now, just use some hard-coded values.
         // (1) Create a pulse library (pulse name to samples)
         PulseLib testPulseLib;
         testPulseLib.emplace("pulse1", std::vector<std::complex<double>> { 0.1, 0.2, 0.1, 0.0, -0.1, -0.2, 0.1, 0.1, 0.05 });
         
         // (2) Backend configs:
         BackendChannelConfigs backendConfig;
         {
            backendConfig.dt = 1.0;
            backendConfig.nb_dChannels = 1;
            backendConfig.loFregs_dChannels = { nu };
            backendConfig.pulseLib = testPulseLib;
         }

         // Pulse schedule entries
         PulseScheduleEntry testPulseScheduleEntry;
         {
            testPulseScheduleEntry.name = "pulse1";
            testPulseScheduleEntry.startTime = 0.0;
            testPulseScheduleEntry.stopTime = stopTime;
         }
         
         PulseScheduleRegistry testPulseSchedule;
         const size_t channelId = 0;
         testPulseSchedule.emplace(channelId, std::vector<PulseScheduleEntry> { testPulseScheduleEntry });
         
         // Initialize the controller
         m_pulseChannelController->Initialize(backendConfig, testPulseSchedule, {});
      }

      XACC_QuaC_InitializePulseSim(buffer->size(), dt, stopTime, stepMax, reinterpret_cast<PulseChannelProvider*>(m_pulseChannelController.get()));
      
      // Debug:
      XACC_QuaC_SetLogVerbosity(DEBUG_DIAG);
      
      // TODO: Parse the QObj to get the Hamiltonian params
      // For now, just hard-coded
      // H = -pi*nu*sigma_z + D(t)*sigma_x
      // D(t) is a Gaussian pulse 
      
      // Time-independent terms:
      XACC_QuaC_AddConstHamiltonianTerm1("Z", 0, { -omega/2.0, 0.0});

      // Time-dependent term:
      XACC_QuaC_AddTimeDependentHamiltonianTerm1("X", 0, "D0");

      XACC_QuaC_AddQubitDecay(0, kappa);
   }

   void PulseVisitor::solve() 
   {
      std::cout << "Pulse simulator: solving the Hamiltonian. \n";
      double* results = nullptr;
      
      TSData* tsData;  
      int nbSteps;
      const auto resultSize = XACC_QuaC_RunPulseSim(&results, &nbSteps, &tsData);
      
#ifdef EXPORT_TS_DATA_AS_CSV
      writeTimesteppingDataToCsv("output", tsData, nbSteps);
#endif

      std::cout << "Final result: ";
      for (int i = 0; i < resultSize; ++i)
      {
         std::cout << results[i] << ", ";
      }
      std::cout << "\n";

      free(results);
   }

   void PulseVisitor::finalize() 
   {     
      std::cout << "Pulse simulator: Finalized. \n";
      XACC_QuaC_Finalize();
   }
}