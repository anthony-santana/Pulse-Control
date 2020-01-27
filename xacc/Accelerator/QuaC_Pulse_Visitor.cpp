#include "QuaC_Pulse_Visitor.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include <chrono>
#include <iomanip>
#include <cassert>
#include "Pulse.hpp"
#include "Scheduler.hpp"

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
      if (in_channelName.front() == 'd' || in_channelName.front() == 'D')
      {
         const std::string channelNum = in_channelName.substr(1);
         return std::make_pair(ChannelType::D, std::stoi(channelNum)); 
      }
      if (in_channelName.front() == 'u' || in_channelName.front() == 'U')
      {
         const std::string channelNum = in_channelName.substr(1);
         return std::make_pair(ChannelType::U, std::stoi(channelNum)); 
      }

      return std::make_pair(ChannelType::Invalid, 0);
   }
   
   double LoadFromPulse(FrameChangeCommandEntry& io_fcCommand, xacc::quantum::Pulse& in_pulseInst, double in_sampleDt)
   {
      assert(in_pulseInst.name() == "fc");
      io_fcCommand.startTime = in_pulseInst.start() * in_sampleDt;
      io_fcCommand.phase = in_pulseInst.getParameter(0).as<double>();
      return io_fcCommand.startTime;
   }

   double LoadFromPulse(PulseScheduleEntry& io_pulseEntry, xacc::quantum::Pulse& in_pulseInst, double in_sampleDt)
   {
      io_pulseEntry.name = in_pulseInst.name();
      io_pulseEntry.startTime = in_pulseInst.start() * in_sampleDt;
      const double expectedStopTime = (in_pulseInst.duration() + in_pulseInst.start()) * in_sampleDt;
      io_pulseEntry.stopTime = expectedStopTime;
      return io_pulseEntry.stopTime;
   }

   // Helper to insert Pulse/FC entries to pulse constroller scheduling map
   template<class PulseItemType>
   void UpdatePulseScheduleMap(std::unordered_map<size_t, std::vector<PulseItemType>>& io_mapToUpdate, xacc::quantum::Pulse& in_pulseInst, PulseChannelController& io_pulseController, double& out_endTime)
   {
      const std::string channelName = in_pulseInst.channel();
      PulseItemType entry;
      out_endTime = LoadFromPulse(entry, in_pulseInst, io_pulseController.GetBackendConfigs().dt);
      const auto channelTypeIdPair = PulseChannelFromString(channelName);
      switch (channelTypeIdPair.first)
      {
         case ChannelType::D:
         {
            const auto channelId = io_pulseController.GetDriveChannelId(channelTypeIdPair.second);
            auto existingScheduleIter = io_mapToUpdate.find(channelId);
            if (existingScheduleIter == io_mapToUpdate.end())
            {
               io_mapToUpdate.emplace(channelId, std::vector<PulseItemType> { entry });
            }
            else
            {
               existingScheduleIter->second.emplace_back(entry);
            }
            break;
         }
         case ChannelType::U:
         {
            const auto channelId = io_pulseController.GetControlChannelId(channelTypeIdPair.second);
            auto existingScheduleIter = io_mapToUpdate.find(channelId);
            if (existingScheduleIter == io_mapToUpdate.end())
            {
               io_mapToUpdate.emplace(channelId, std::vector<PulseItemType> { entry });
            }
            else
            {
               existingScheduleIter->second.emplace_back(entry);
            }
            break;
         } 
         case ChannelType::Invalid:
            xacc::error("Invalid Pulse channel named '" + channelName + "'\n");
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
      double kappa = 1e-64;
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
         // XACC_QuaC_SetLogVerbosity(DEBUG_DIAG);
      }

      {
         // Step 2: set up the Hamiltonian
         if (in_params.stringExists("hamiltonian"))
         {
            const auto hamiltonianJson = in_params.getString("hamiltonian");         
            auto parser = xacc::getService<QuaC::HamiltonianParsingUtil>("default");
            
            std::function<void(QuaC::HamiltonianTerm&)> iterFn = [&](QuaC::HamiltonianTerm& in_term) -> void {
               in_term.apply(this);
            };

            if (!parser->tryParse(hamiltonianJson, iterFn))
            {
               xacc::error("Failed to parse the Hamiltonian!");
            }

            for (int i = 0; i < buffer->size(); ++i)
            {
               XACC_QuaC_AddQubitDecay(i, kappa);
            }
         }
         else
         {
           // TODO: clean up this code
           // Some special case handling, to be removed
            if (backendName == "Fake1Q")
            {
               // The Hamiltonian is (use the one from arXiv): 
               // H = -pi*v0*sigma_z + D(t)*sigma_x  
               // Time-independent terms:
               assert(backendConfig.loFregs_dChannels.size() == 1);
               XACC_QuaC_AddConstHamiltonianTerm1("Z", 0, { -M_PI * backendConfig.loFregs_dChannels[0], 0.0});    
               
               // Time-dependent term (drive channel 0):
               // Note: we calibrate this fake one qubit system to work well with the QObj data that we have
               // from another device. 
               const double tdCoeff = 1.5; 
               XACC_QuaC_AddTimeDependentHamiltonianTerm1("X", 0, m_pulseChannelController->GetDriveChannelId(0), tdCoeff);
               
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
                  // Reference: https://github.com/Qiskit/qiskit-terra/blob/13bc243364553667f6410b9a2f7a315c90bb598f/qiskit/test/mock/fake_openpulse_2q.py
                  const double r = 0.02;
                  // Drive channel D0: 2*np.pi*r*X0||D0
                  XACC_QuaC_AddTimeDependentHamiltonianTerm1("X", 0, m_pulseChannelController->GetDriveChannelId(0), 2 * M_PI * r);
                  
                  // D1: 2*np.pi*r*X1||D1
                  XACC_QuaC_AddTimeDependentHamiltonianTerm1("X", 1, m_pulseChannelController->GetDriveChannelId(1), 2 * M_PI * r);
                  
                  // U0: 2*np.pi*r*X1||U0
                  XACC_QuaC_AddTimeDependentHamiltonianTerm1("X", 1, m_pulseChannelController->GetControlChannelId(0), 2 * M_PI * r);

                  // U1: 2*np.pi*r*X0||U1              
                  XACC_QuaC_AddTimeDependentHamiltonianTerm1("X", 0, m_pulseChannelController->GetControlChannelId(1), 2 * M_PI * r);
               }

               {
                  // Add some decay
                  XACC_QuaC_AddQubitDecay(0, kappa);     
                  XACC_QuaC_AddQubitDecay(1, kappa);     
               }
            }        
         }         
      }
   }

   int PulseVisitor::GetChannelId(const std::string& in_channelName) 
   {
      const auto channelTypeIdPair = PulseChannelFromString(in_channelName);
      switch (channelTypeIdPair.first)
      {
         case ChannelType::D: return m_pulseChannelController->GetDriveChannelId(channelTypeIdPair.second);
         case ChannelType::U: return m_pulseChannelController->GetControlChannelId(channelTypeIdPair.second);
      }

      return -1;
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

   void PulseVisitor::schedulePulses(const std::shared_ptr<CompositeInstruction>& in_pulseInstruction) 
   {
      // Check if any raw pulse instructions are referring to pulse from the library,
      // i.e. not carrying correct duration (sample size) data
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
            
            const std::string pulseName = pulse->name();
            // This pulse doesn't have duration data 
            // i.e. perhaps just refers to a known pulse, hence update its data accordingly.
            if (pulseName != "fc" && pulseName != "acquire" && pulse->duration() == 0)
            {
               const size_t pulseSampleSizeFromLib = m_pulseChannelController->GetBackendConfigs().getPulseSampleSize(pulseName);
               if (pulseSampleSizeFromLib > 0)
               {
                  // Update the pulse data
                  pulse->setDuration(pulseSampleSizeFromLib);
               }
            }
         }
      }

      // After we have adjust the pulse duration (if missing)
      // Use XACC Pulse scheduler to schedule pulses:
      // Debug:
      // std::cout << "Before Scheduled : \n" << in_pulseInstruction->toString() << "\n";
   
      auto scheduler = xacc::getService<xacc::Scheduler>("pulse");
      scheduler->schedule(in_pulseInstruction);
      
      // Debug:
      // std::cout << "After Scheduled : \n" << in_pulseInstruction->toString() << "\n";
   }

   void PulseVisitor::solve(const std::shared_ptr<CompositeInstruction>& in_pulseInstruction) 
   {
      // Step 1: schedule the pulse program
      schedulePulses(in_pulseInstruction);      
      
      // Step 2: initialize the pulse channel controller with the scheduled pulses
      PulseScheduleRegistry allPulseSchedules;
      FrameChangeScheduleRegistry allFcSchedules;
      double simStopTime = 0.0;
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

            // Handle Frame changes
            if (pulseName == "fc")
            {
               double endTime = 0.0;
               UpdatePulseScheduleMap(allFcSchedules, *pulse, *m_pulseChannelController, endTime);
               if (endTime > simStopTime)
               {
                  simStopTime = endTime;
               }
            }
            else if (pulseName == "acquire")
            {
               // TODO: we don't support acquire yet
               xacc::warning("Acquire is not supported. Skipped!!!!");
            }
            else
            {
               // This is a pulse instruction (pulse name + sample)
               
               // Check in case this is a new pulse that may be constructed on-the-fly
               // e.g. we can, at IR-level, create a new pulse (arbitrary name and supply the samples),
               // then add to the composite. We will gracefully handle that by adding that user-defined pulse to the library. 
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
               
               double endTime = 0.0;
               UpdatePulseScheduleMap(allPulseSchedules, *pulse, *m_pulseChannelController, endTime);
               if (endTime > simStopTime)
               {
                  simStopTime = endTime;
               }
            }
         }
      }
      
      // Initialize the channel controller with pulse and fc commands     
      m_pulseChannelController->Initialize(allPulseSchedules, allFcSchedules);

      std::cout << "Pulse simulator: solving the Hamiltonian. \n";
      double* results = nullptr;
      
      // Note: This dt is the solver step size (which may be adaptive),
      // this should be smaller than the Pulse dt (sample step size).
      double dt = m_pulseChannelController->GetBackendConfigs().dt/1000.0;
      
      // Add some extra time (device dt) to simulation time
      simStopTime += m_pulseChannelController->GetBackendConfigs().dt; 
    
      int stepMax = static_cast<int>(2.0 * std::ceil(simStopTime/dt));
      
      TSData* tsData;  
      int nbSteps;
      
      // Step 3: Run the simulation
      const auto resultSize = XACC_QuaC_RunPulseSim(dt, simStopTime, stepMax, &results, &nbSteps, &tsData);
      
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