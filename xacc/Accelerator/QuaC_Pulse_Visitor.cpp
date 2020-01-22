#include "QuaC_Pulse_Visitor.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include <chrono>
#include <iomanip>
#include <cassert>
#include "Pulse.hpp"

extern "C" {
#include "interface_xacc_ir.h"
}

#ifdef ENABLE_MOCKING
#include "FakeOpenPulse1Q.hpp"
#include "FakeOpenPulse2Q.hpp"
#endif

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

      const auto nbChannels = in_tsData[0].nbChannels;      
      const auto nbPopulations = in_tsData[0].nbPops;
      // Header
      outputFile << "Time, "; 
      for(int j = 0; j < nbChannels; ++j)
      {
         outputFile << "Channel[" << j << "], ";
      }     
      for(int j = 0; j < nbPopulations; ++j)
      {
         outputFile << "Population[" << j << "], ";
      }
      outputFile << "\n";

      // Data
      for (int i = 0; i < in_nbSteps; ++i)
      {
         const TSData dataAtStep = in_tsData[i];
         // First column is time:
         outputFile << dataAtStep.time << ", ";
         // Channel data columns
         for(int j = 0; j < nbChannels; ++j)
         {
            outputFile << dataAtStep.channelData[j] << ", ";
         }
         // Population columns
         for(int j = 0; j < nbPopulations; ++j)
         {
            outputFile << dataAtStep.populations[j] << ", ";
         }
         outputFile << "\n";
      }
      outputFile.close();
      std::cout << "Time-stepping data is written to file '" << fileName << "'\n";
   } 

   enum class ChannelType { Invalid, D, U };
   std::pair<ChannelType, int> PulseChannelFromString(const std::string& in_channelName)
   {
      if (in_channelName.front() != 'd' && in_channelName.front() != 'u')
      {
         return std::make_pair(ChannelType::Invalid, 0);
      }

      if (in_channelName.front() == 'd')
      {
         const std::string channelNum = in_channelName.substr(1);
         return std::make_pair(ChannelType::D, std::stoi(channelNum)); 
      }
      if (in_channelName.front() == 'u')
      {
         const std::string channelNum = in_channelName.substr(1);
         return std::make_pair(ChannelType::U, std::stoi(channelNum)); 
      }
   }
}

// Export the time-stepping data to csv (e.g. for plotting)
// Uncomment to get the data exported
#define EXPORT_TS_DATA_AS_CSV

namespace QuaC {   
   void PulseVisitor::initialize(std::shared_ptr<AcceleratorBuffer> buffer, const HeterogeneousMap& in_params, const PulseLib& in_importedPulses) 
   {
      // Debug
      std::cout << "Initialize Pulse simulator \n";
         
      // Qubit decay: just use a very small value
      // TODO: we can convert the T1 data from backend data to this param
      double kappa = 0.0001;
     
     
      double nu = 5.0;
      double omega = 2 * M_PI * nu;

      BackendChannelConfigs backendConfig;
      std::string backendName;
      if(in_params.stringExists("backend"))
      {
         backendName = in_params.getString("backend");         
         if (backendName == "Fake1Q")
         {
            FakePulse1Q fakePulse;
            backendConfig = fakePulse.backendConfig;
         }
         else if (backendName == "Fake2Q")
         {
            FakePulse2Q fakePulse;
            backendConfig = fakePulse.backendConfig;
         }
         else
         {
            xacc::error("Unknown backend named '" + backendName + "'.\n");
         }
      } 
      else 
      {
         xacc::error("No backend specified. Exiting...\n");
      }
      
      // ===================================================================================================   
      // Test code: Import those pulses from ibmq_poughkeepsie backend for testing.
      // Note: we just import the pulse library (i.e. list of samples), not simulate its whole Hamiltonian.
      //====================================================================================================
      // We should eventually import everything from Json file, for now, we mock the dynamic (i.e. Hamiltonian)
      // but import real pulse library to test the integration.
      for (const auto& importedPulse : in_importedPulses)
      {
         backendConfig.addOrReplacePulse(importedPulse.first, importedPulse.second);
      }
      //===============================================================================================


      // Initialize the pulse constroller:
      // Create a pulse controller
      m_pulseChannelController = std::make_unique<PulseChannelController>(backendConfig);
      
      {
         // Step 1: Initialize the QuaC solver
         XACC_QuaC_InitializePulseSim(buffer->size(), reinterpret_cast<PulseChannelProvider*>(m_pulseChannelController.get()));
         // Debug:
         XACC_QuaC_SetLogVerbosity(DEBUG_DIAG);
      }

      {
         // Step 2: set up the Hamiltonian
         // TODO: Parse the QObj to get the Hamiltonian params
         // For now, just do it manually for the 1 and 2 Q hamiltonian
         if (backendName == "Fake1Q")
         {
            // The Hamiltonian is (use the one from arXiv): 
            // H = -pi*v0*sigma_z + D(t)*sigma_x  
            // Time-independent terms:
            assert(backendConfig.loFregs_dChannels.size() == 1);
            XACC_QuaC_AddConstHamiltonianTerm1("Z", 0, { -M_PI * backendConfig.loFregs_dChannels[0], 0.0});    
            
            // Time-dependent term (drive channel 0):
            XACC_QuaC_AddTimeDependentHamiltonianTerm1("X", 0, m_pulseChannelController->GetDriveChannelId(0));
            
            // Add some decay
            XACC_QuaC_AddQubitDecay(0, kappa);     
         }
         else if (backendName == "Fake2Q")
         {
            // The Hamiltonian is: 
            // Notation: 
            // {'X': sigmax, 'Y': sigmay, 'Z': sigmaz,
            //   'Sp': create, 'Sm': destroy, 'I': qeye,
            //   'O': num, 'P': project, 'A': destroy,
            //   'C': create, 'N': num}
            // ["np.pi*(2*v0-alpha0)*O0", "np.pi*alpha0*O0*O0", "2*np.pi*r*X0||D0",
            //  "2*np.pi*r*X0||U1", "2*np.pi*r*X1||U0", "np.pi*(2*v1-alpha1)*O1",
            //  "np.pi*alpha1*O1*O1", "2*np.pi*r*X1||D1", "2*np.pi*j*(Sp0*Sm1+Sm0*Sp1)"],
            
            // Time-independent terms:
            {
               const double v0 = 5.00;
               const double v1 = 5.10;
               const double alpha0 = -0.33;
               const double alpha1 = -0.33;
               const double j = 0.01;
               XACC_QuaC_AddConstHamiltonianTerm1("O", 0, { M_PI * (2*v0-alpha0), 0.0 });
               XACC_QuaC_AddConstHamiltonianTerm2("O", 0, "O", 0, { M_PI * alpha0, 0.0});        
               XACC_QuaC_AddConstHamiltonianTerm1("O", 1, { M_PI * (2*v1-alpha1), 0.0 });
               XACC_QuaC_AddConstHamiltonianTerm2("O", 1, "O", 1, { M_PI * alpha1, 0.0});      
               XACC_QuaC_AddConstHamiltonianTerm2("SP", 0, "SM", 1, { 2.0 * M_PI * j, 0.0});      
               XACC_QuaC_AddConstHamiltonianTerm2("SM", 0, "SP", 1, { 2.0 * M_PI * j, 0.0});
            }
            // Time-dependent terms
            // TODO: we should add the multiplication coefficient to the API 
            {
               // Drive channel D0: 2*np.pi*r*X0||D0
               XACC_QuaC_AddTimeDependentHamiltonianTerm1("X", 0, m_pulseChannelController->GetDriveChannelId(0));
               
               // D1: 2*np.pi*r*X1||D1
               XACC_QuaC_AddTimeDependentHamiltonianTerm1("X", 1, m_pulseChannelController->GetDriveChannelId(1));
               
               // U0: 2*np.pi*r*X1||U0
               XACC_QuaC_AddTimeDependentHamiltonianTerm1("X", 1, m_pulseChannelController->GetControlChannelId(0));

               // U1: 2*np.pi*r*X0||U1              
               XACC_QuaC_AddTimeDependentHamiltonianTerm1("X", 0, m_pulseChannelController->GetControlChannelId(1));
            }

            {
               // Add some decay
               XACC_QuaC_AddQubitDecay(0, kappa);     
               XACC_QuaC_AddQubitDecay(1, kappa);     
            }
         }        
      }
   }


   std::vector<std::complex<double>> PulseVisitor::PulseSamplesToComplexVec(const std::vector<std::vector<double>>& in_samples)
   {
      std::vector<std::complex<double>> pulseSamples;
      pulseSamples.reserve(in_samples.size());
      for (const auto& sample : in_samples)
      {
         if (sample.empty())
         {
               pulseSamples.emplace_back(0.0);
         }
         else if (sample.size() == 1)
         {
               pulseSamples.emplace_back(sample[0]);
         }
         else if (sample.size() == 2)
         {
               pulseSamples.emplace_back(std::complex<double>{ sample[0], sample[1] });
         }
         else
         {
            // Malformed
            std::cout << "Invalid data encountered while processing pulse samples.\n";
            return {};
         }                        
      }

      return pulseSamples;
   }


   void PulseVisitor::solve(const std::shared_ptr<CompositeInstruction>& in_pulseInstruction) 
   {
      // Add *compiled* pulse instructions to the pulse controller's schedule
      // i.e. the pulse controller will drive all time-dependent channels according to the pulse CompositeInstruction
      std::unordered_map<std::string, double> channelToTime;
      PulseScheduleRegistry allPulseSchedules;

      xacc::InstructionIterator it(in_pulseInstruction);
      while (it.hasNext()) 
      {
         auto nextInst = it.next();
         if (nextInst->isEnabled() && !nextInst->isComposite()) 
         {
            auto pulse = std::dynamic_pointer_cast<xacc::quantum::Pulse>(nextInst);
            if (!pulse) 
            {
               xacc::error("Invalid instruction in pulse program.");
            }
            const std::string channelName = pulse->channel();
            const std::string pulseName = pulse->name();
            if (!m_pulseChannelController->GetBackendConfigs().hasPulseName(pulseName))
            {
               // Perhaps it's a dynamic pulse (constructed on the fly using IR-level API)
               if (pulse->getSamples().empty())
               {
                  xacc::error("Invalid instruction in pulse program.");
               }
               // Imported the pulse (on-demand)
               m_pulseChannelController->GetBackendConfigs().addOrReplacePulse(pulseName, PulseSamplesToComplexVec(pulse->getSamples()));
            }  

            if (channelToTime.find(channelName) == channelToTime.end())
            {
               channelToTime.emplace(channelName, 0);
            }    

            auto& currentChannelTime = channelToTime[channelName];
            PulseScheduleEntry scheduleEntry;
            scheduleEntry.name = pulseName;
            scheduleEntry.startTime = currentChannelTime;
            const double expectedStopTime = m_pulseChannelController->GetBackendConfigs().getPulseDuration(pulseName) + currentChannelTime;
            scheduleEntry.stopTime = expectedStopTime;
            currentChannelTime = expectedStopTime;

            const auto channelTypeIdPair = PulseChannelFromString(channelName);
            switch (channelTypeIdPair.first)
            {
               case ChannelType::D:
               {
                  const auto channelId = m_pulseChannelController->GetDriveChannelId(channelTypeIdPair.second);
                  auto existingScheduleIter = allPulseSchedules.find(channelId);
                  if (existingScheduleIter == allPulseSchedules.end())
                  {
                     allPulseSchedules.emplace(channelId, std::vector<PulseScheduleEntry> { scheduleEntry });
                  }
                  else
                  {
                     existingScheduleIter->second.emplace_back(scheduleEntry);
                  }
                  break;
              }
              case ChannelType::U:
              {
                  const auto channelId = m_pulseChannelController->GetControlChannelId(channelTypeIdPair.second);
                  auto existingScheduleIter = allPulseSchedules.find(channelId);
                  if (existingScheduleIter == allPulseSchedules.end())
                  {
                     allPulseSchedules.emplace(channelId, std::vector<PulseScheduleEntry> { scheduleEntry });
                  }
                  else
                  {
                     existingScheduleIter->second.emplace_back(scheduleEntry);
                  }
                  break;
              } 
              case ChannelType::Invalid:
                  xacc::error("Invalid Pulse channel named '" + channelName + "'\n");
            }
         }
      }
      
      
      // TODO: Implement phase changes and other commands      
      m_pulseChannelController->Initialize(allPulseSchedules, {});

      std::cout << "Pulse simulator: solving the Hamiltonian. \n";
      double* results = nullptr;
      
      // Note: This dt is the solver step size (which may be adaptive),
      // this should be smaller than the Pulse dt (sample step size).
      double dt = m_pulseChannelController->GetBackendConfigs().dt/100.0;
      double stopTime = 0.0;
      for (const auto& kv : channelToTime) 
      {
         if (kv.second > stopTime)
         {
            // Stop time is the max time on all channels
            stopTime = kv.second;
         }
      } 
      // Add some extra time 
      stopTime += 1.0; 
    
      int stepMax = static_cast<int>(std::ceil(stopTime/dt));
      
      TSData* tsData;  
      int nbSteps;
      const auto resultSize = XACC_QuaC_RunPulseSim(dt, stopTime, stepMax, &results, &nbSteps, &tsData);
      
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