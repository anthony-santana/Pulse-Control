#include "QuaC_Pulse_Visitor.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include <chrono>
#include <iomanip>
#include <cassert>
#include "Scheduler.hpp"
#include "PulseSystemModel.hpp"

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

   std::string constructPulseCommandDef(xacc::quantum::Gate& in_gate)
   {
      const auto getGateCommandDefName = [](xacc::quantum::Gate& in_gate) -> std::string {
         static std::locale loc;
         const std::string gateName = std::tolower(in_gate.name(), loc);
         if (gateName == "cnot")
         {
            return "cx";
         }
         else
         {
            return gateName;
         }
      };

      std::string result = "pulse::" + getGateCommandDefName(in_gate);

      for (const auto& qIdx : in_gate.bits())
      {
         result += ("_" + std::to_string(qIdx));
      }

      return result;      
   }

   std::shared_ptr<xacc::CompositeInstruction> constructU3CmdDefComposite(size_t qIdx)
   {
      const auto cmdDefName = "pulse::u3_" + std::to_string(qIdx);
      return xacc::ir::asComposite(xacc::getContributedService<xacc::Instruction>(cmdDefName));
   }
}

// Export the time-stepping data to csv (e.g. for plotting)
// Uncomment to get the data exported
#define EXPORT_TS_DATA_AS_CSV

namespace QuaC {   
   void PulseVisitor::initialize(std::shared_ptr<AcceleratorBuffer> buffer, PulseSystemModel* in_systemModel, const HeterogeneousMap& in_params) 
   {
      // Debug
      std::cout << "Initialize Pulse simulator \n";
      auto provider = xacc::getIRProvider("quantum");
      // Should create a hash name here:
      m_pulseComposite = provider->createComposite("PulseComposite");
      m_systemModel = in_systemModel;

      // Initialize the pulse constroller:
      // Create a pulse controller
      m_pulseChannelController = std::make_unique<PulseChannelController>(m_systemModel->getChannelConfigs());
      
      {
         // Step 1: Initialize the QuaC solver
         XACC_QuaC_InitializePulseSim(buffer->size(), reinterpret_cast<PulseChannelProvider*>(m_pulseChannelController.get()));
         // Debug:
         // XACC_QuaC_SetLogVerbosity(DEBUG_DIAG);
      }

      {
         // Step 2: set up the Hamiltonian
         for (const auto& term : m_systemModel->getHamiltonian().getTerms())
         {
            term->apply(this);
         }

         for (int i = 0; i < buffer->size(); ++i)
         {
            // TODO: try to track down why QuaC is throwing exceptions if there is no decay!!!
            // For now, if none decay (T1) specified, just use a super small value.
            const double DEFAULT_KAPPA = 1e-64;
            const double kappa = m_systemModel->getQubitT1(i) == 0.0 ? DEFAULT_KAPPA : (1.0 / m_systemModel->getQubitT1(i));
            XACC_QuaC_AddQubitDecay(i, kappa);
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
      std::cout << "Before Scheduled : \n" << in_pulseInstruction->toString() << "\n";
   
      auto scheduler = xacc::getService<xacc::Scheduler>("pulse");
      scheduler->schedule(in_pulseInstruction);
      
      // Debug:
      std::cout << "After Scheduled : \n" << in_pulseInstruction->toString() << "\n";
   }

   void PulseVisitor::solve() 
   {
      // Step 1: schedule the pulse program
      schedulePulses(m_pulseComposite);      
      
      // Step 2: initialize the pulse channel controller with the scheduled pulses
      PulseScheduleRegistry allPulseSchedules;
      FrameChangeScheduleRegistry allFcSchedules;
      double simStopTime = 0.0;
      xacc::InstructionIterator it(m_pulseComposite);
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

   void PulseVisitor::visit(Hadamard& h)  
   {
      const auto commandDef = constructPulseCommandDef(h);

      if (xacc::hasContributedService<xacc::Instruction>(commandDef))
      {
         auto pulseInst = xacc::getContributedService<xacc::Instruction>(commandDef);
         m_pulseComposite->addInstruction(pulseInst);
      }
      else
      {
         // H = U2(0, pi) = U3(pi/2, 0, pi)
         const auto asU3 = constructU3CmdDefComposite(h.bits()[0]);
         auto hCmdDef = (*asU3)({xacc::constants::pi/2.0, 0.0, xacc::constants::pi});
         m_pulseComposite->addInstruction(hCmdDef);
      }
   }
   
   void PulseVisitor::visit(CNOT& cnot)  
   {
      const auto commandDef = constructPulseCommandDef(cnot);

      if (xacc::hasContributedService<xacc::Instruction>(commandDef))
      {
         auto pulseInst = xacc::getContributedService<xacc::Instruction>(commandDef);
         m_pulseComposite->addInstruction(pulseInst);
      }
      else
      {
        xacc::error("CNOT pulse cmd-def is not provided!");
      }
   }
   
   void PulseVisitor::visit(Rz& rz)  
   {
      const auto commandDef = constructPulseCommandDef(rz);

      if (xacc::hasContributedService<xacc::Instruction>(commandDef))
      {
         auto pulseInst = xacc::getContributedService<xacc::Instruction>(commandDef);
         m_pulseComposite->addInstruction(pulseInst);
      }
      else
      {
         // Rz(theta) = U1(theta) = U3(0, 0, theta)
         const auto asU3 = constructU3CmdDefComposite(rz.bits()[0]);
         auto rzCmdDef = (*asU3)({0.0, 0.0, rz.getParameter(0).as<double>()});
         m_pulseComposite->addInstruction(rzCmdDef);
      }
   }
   
   void PulseVisitor::visit(Ry& ry)  
   {
      const auto commandDef = constructPulseCommandDef(ry);

      if (xacc::hasContributedService<xacc::Instruction>(commandDef))
      {
         auto pulseInst = xacc::getContributedService<xacc::Instruction>(commandDef);
         m_pulseComposite->addInstruction(pulseInst);
      }
      else
      {
         // Ry(theta) = U3(theta, 0, 0)
         const auto asU3 = constructU3CmdDefComposite(ry.bits()[0]);
         auto ryCmdDef = (*asU3)({ry.getParameter(0).as<double>(), 0.0, 0.0});
         m_pulseComposite->addInstruction(ryCmdDef);
      }
   }
   
   void PulseVisitor::visit(Rx& rx)  
   {
      const auto commandDef = constructPulseCommandDef(rx);

      if (xacc::hasContributedService<xacc::Instruction>(commandDef))
      {
         auto pulseInst = xacc::getContributedService<xacc::Instruction>(commandDef);
         m_pulseComposite->addInstruction(pulseInst);
      }
      else
      {
         // Rx(theta) = U3(theta, -pi/2, pi/2)
         const auto asU3 = constructU3CmdDefComposite(rx.bits()[0]);
         auto rxCmdDef = (*asU3)({rx.getParameter(0).as<double>(), -xacc::constants::pi/2.0, xacc::constants::pi/2.0});
         m_pulseComposite->addInstruction(rxCmdDef);
      }
   }
   
   void PulseVisitor::visit(X& x)  
   {
      const auto commandDef = constructPulseCommandDef(x);

      if (xacc::hasContributedService<xacc::Instruction>(commandDef))
      {
         auto pulseInst = xacc::getContributedService<xacc::Instruction>(commandDef);
         m_pulseComposite->addInstruction(pulseInst);
      }
      else
      {
         // X = U3(pi, 0, pi)
         const auto asU3 = constructU3CmdDefComposite(x.bits()[0]);
         auto xCmdDef = (*asU3)({xacc::constants::pi, 0.0, xacc::constants::pi});
         m_pulseComposite->addInstruction(xCmdDef);
      }
   }
   
   void PulseVisitor::visit(Y& y)  
   {
      const auto commandDef = constructPulseCommandDef(y);

      if (xacc::hasContributedService<xacc::Instruction>(commandDef))
      {
         auto pulseInst = xacc::getContributedService<xacc::Instruction>(commandDef);
         m_pulseComposite->addInstruction(pulseInst);
      }
      else
      {
         // Y = U3(pi, pi/2, pi/2)
         const auto asU3 = constructU3CmdDefComposite(y.bits()[0]);
         auto yCmdDef = (*asU3)({xacc::constants::pi, xacc::constants::pi/2.0, xacc::constants::pi/2.0});
         m_pulseComposite->addInstruction(yCmdDef);
      }
   }
   
   void PulseVisitor::visit(Z& z)  
   {
      const auto commandDef = constructPulseCommandDef(z);

      if (xacc::hasContributedService<xacc::Instruction>(commandDef))
      {
         auto pulseInst = xacc::getContributedService<xacc::Instruction>(commandDef);
         m_pulseComposite->addInstruction(pulseInst);
      }
      else
      {
         // Z = U1(pi) = U3(0, 0, pi)
         const auto asU3 = constructU3CmdDefComposite(z.bits()[0]);
         auto zCmdDef = (*asU3)({0.0, 0.0, xacc::constants::pi});
         m_pulseComposite->addInstruction(zCmdDef);
      }
   }
   
   void PulseVisitor::visit(CY& cy)  
   {
      // CY(a,b) = sdg(b); cx(a,b); s(b); 
      auto sdg = std::make_shared<Sdg>(cy.bits()[1]);
      auto cx = std::make_shared<CNOT>(cy.bits());
      auto s = std::make_shared<S>(cy.bits()[1]);
      visit(*sdg);
      visit(*cx);
      visit(*s);
   }
   
   void PulseVisitor::visit(CZ& cz)  
   {
      // CZ(a,b) =  H(b); CX(a,b); H(b); 
      auto h1 = std::make_shared<Hadamard>(cz.bits()[1]);
      auto cx = std::make_shared<CNOT>(cz.bits());
      auto h2 = std::make_shared<Hadamard>(cz.bits()[0]);
      visit(*h1);
      visit(*cx);
      visit(*h2);
   }
   
   void PulseVisitor::visit(Swap& s)  
   {
     // SWAP(a,b) =  CX(a,b); CX(b,a); CX(a,b); 
      auto cx1 = std::make_shared<CNOT>(s.bits());
      auto cx2 = std::make_shared<CNOT>(s.bits()[1], s.bits()[0]);
      auto cx3 = std::make_shared<CNOT>(s.bits());
      visit(*cx1);
      visit(*cx2);
      visit(*cx3);
   }
   
   void PulseVisitor::visit(CRZ& crz)  
   {
      // CRZ(theta)(a,b) = U3(0, 0, theta/2)(b); CX(a,b); U3(0, 0, -theta/2)(b); CX(a,b);
      const double theta = crz.getParameter(0).as<double>();
      {
         const auto asU3 = constructU3CmdDefComposite(crz.bits()[1]);
         auto cmdDef = (*asU3)({0.0, 0.0, theta/2.0});
         m_pulseComposite->addInstruction(cmdDef);
      }
      {
         auto cx = std::make_shared<CNOT>(crz.bits());
         visit(*cx);
      }
      {
         const auto asU3 = constructU3CmdDefComposite(crz.bits()[1]);
         auto cmdDef = (*asU3)({0.0, 0.0, -theta/2.0});
         m_pulseComposite->addInstruction(cmdDef);
      }
      {
         auto cx = std::make_shared<CNOT>(crz.bits());
         visit(*cx);
      }
   }
   
   void PulseVisitor::visit(CH& ch)  
   {
      // CH(a,b) = S(b); H(b); T(b); CX(a,b); Tdg(b); H(b); Sdg(b); 
      {
         auto s = std::make_shared<S>(ch.bits()[1]);
         visit(*s);
      }
      {
         auto h = std::make_shared<Hadamard>(ch.bits()[1]);
         visit(*h);
      }
      {
         auto t = std::make_shared<T>(ch.bits()[1]);
         visit(*t);
      }
      {
         auto cx = std::make_shared<CNOT>(ch.bits());
         visit(*cx);
      }
      {
         auto tdg = std::make_shared<Tdg>(ch.bits()[1]);
         visit(*tdg);
      }
      {
         auto h = std::make_shared<Hadamard>(ch.bits()[1]);
         visit(*h);
      }
      {
         auto sdg = std::make_shared<Sdg>(ch.bits()[1]);
         visit(*sdg);
      }
   }
   
   void PulseVisitor::visit(S& s)  
   {
      const auto commandDef = constructPulseCommandDef(s);

      if (xacc::hasContributedService<xacc::Instruction>(commandDef))
      {
         auto pulseInst = xacc::getContributedService<xacc::Instruction>(commandDef);
         m_pulseComposite->addInstruction(pulseInst);
      }
      else
      {
         // S = U1(pi/2) = U3(0,0,pi/2)
         const auto asU3 = constructU3CmdDefComposite(s.bits()[0]);
         auto sCmdDef = (*asU3)({0.0, 0.0, xacc::constants::pi/2.0});
         m_pulseComposite->addInstruction(sCmdDef);
      }
   }
   
   void PulseVisitor::visit(Sdg& sdg)  
   {
      const auto commandDef = constructPulseCommandDef(sdg);

      if (xacc::hasContributedService<xacc::Instruction>(commandDef))
      {
         auto pulseInst = xacc::getContributedService<xacc::Instruction>(commandDef);
         m_pulseComposite->addInstruction(pulseInst);
      }
      else
      {
         // S-dagger = U1(-pi/2) = U3(0,0,-pi/2)
         const auto asU3 = constructU3CmdDefComposite(sdg.bits()[0]);
         auto sdgCmdDef = (*asU3)({0.0, 0.0, -xacc::constants::pi/2.0});
         m_pulseComposite->addInstruction(sdgCmdDef);
      }
   }
   
   void PulseVisitor::visit(T& t)  
   {
      const auto commandDef = constructPulseCommandDef(t);

      if (xacc::hasContributedService<xacc::Instruction>(commandDef))
      {
         auto pulseInst = xacc::getContributedService<xacc::Instruction>(commandDef);
         m_pulseComposite->addInstruction(pulseInst);
      }
      else
      {
         // T = U1(pi/4) = U3(0,0,pi/4)
         const auto asU3 = constructU3CmdDefComposite(t.bits()[0]);
         auto tCmdDef = (*asU3)({0.0, 0.0, xacc::constants::pi/4.0});
         m_pulseComposite->addInstruction(tCmdDef);
      }
   }
   
   void PulseVisitor::visit(Tdg& tdg)  
   {
      const auto commandDef = constructPulseCommandDef(tdg);

      if (xacc::hasContributedService<xacc::Instruction>(commandDef))
      {
         auto pulseInst = xacc::getContributedService<xacc::Instruction>(commandDef);
         m_pulseComposite->addInstruction(pulseInst);
      }
      else
      {
         // T-dagger = U1(-pi/4) = U3(0,0,-pi/4)
         const auto asU3 = constructU3CmdDefComposite(tdg.bits()[0]);
         auto tdgCmdDef = (*asU3)({0.0, 0.0, -xacc::constants::pi/4.0});
         m_pulseComposite->addInstruction(tdgCmdDef);
      }
   }
   
   void PulseVisitor::visit(CPhase& cphase)  
   {
      // TODO
      xacc::error("Not supported in Pulse mode.");
   }

   void PulseVisitor::visit(Identity& i)  
   {
      // Ignore for now
      xacc::warning("Ignore Identity gate in Pulse mode.");
   }
   
   void PulseVisitor::visit(U& u)  
   {
      const auto asU3 = constructU3CmdDefComposite(u.bits()[0]);
      // Just pass the params to the cmddef.
      std::vector<double> params;
      for (const auto& param: u.getParameters())
      {
         params.emplace_back(param.as<double>());
      }

      auto uCmdDef = (*asU3)(params);
      m_pulseComposite->addInstruction(uCmdDef);
   }

   void PulseVisitor::visit(Measure& measure)  
   {
      xacc::error("Measure is not supported atm!");
   }

   void PulseVisitor::visit(Pulse& p)  
   {
      // It's already a pulse:
      m_pulseComposite->addInstruction(p.clone());
   }

   void PulseVisitor::finalize() 
   {     
      std::cout << "Pulse simulator: Finalized. \n";
      XACC_QuaC_Finalize();
   }
}