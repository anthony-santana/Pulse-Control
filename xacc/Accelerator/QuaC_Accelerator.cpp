#include "QuaC_Accelerator.hpp"
#include "QuaC_Circuit_Visitor.hpp"
#include "QuaC_Pulse_Visitor.hpp"
#include "json.hpp"
#include "Pulse.hpp"

namespace QuaC {
    void QuaC_Accelerator::initialize(const HeterogeneousMap& params)  
    {
        // DEBUG
        std::cout << ">> DEBUG: QuaC_Accelerator initialize ...\n";
        if (params.stringExists("sim-mode")) 
        {
            const auto requestedMode = params.getString("sim-mode");
            if (requestedMode == "Pulse")
            {
                m_isPulse = true;
                if (params.stringExists("config-json-path"))
                {
                    const auto jsonFileName = params.getString("config-json-path");
                    std::ifstream backendFile(jsonFileName);
                    std::string jjson((std::istreambuf_iterator<char>(backendFile)), std::istreambuf_iterator<char>());
                    // Add those pulses in the library as intructions.
                    contributeInstructions(jjson);
                }
            }
        }

        m_params = params;
        
        if (m_isPulse)
        {
            m_pulseVisitor = std::make_shared<PulseVisitor>();
        }
        else
        {
            m_visitor = std::make_shared<CircuitVisitor>();   
            // TODO
            // Initialize noise, gate time params etc. since QuaC is able to simulate noisy circuits
        }
    }

    void QuaC_Accelerator::execute(std::shared_ptr<AcceleratorBuffer> buffer, const std::shared_ptr<CompositeInstruction> compositeInstruction)  
    {
        if (!m_isPulse)
        {
            m_visitor->initialize(buffer);
            // Walk the IR tree, and visit each node
            InstructionIterator it(compositeInstruction);
            while (it.hasNext()) 
            {
                auto nextInst = it.next();
                if (nextInst->isEnabled()) 
                {
                    nextInst->accept(m_visitor);
                }
            }

            m_visitor->finalize();
        }
        else
        {            
            m_pulseVisitor->initialize(buffer, m_params, m_importedPulses);
            m_pulseVisitor->solve(compositeInstruction);
            m_pulseVisitor->finalize();
        }       
    }
    void QuaC_Accelerator::execute(std::shared_ptr<AcceleratorBuffer> buffer, const std::vector<std::shared_ptr<CompositeInstruction>> compositeInstructions)  
    {
       // TODO
    }

    void QuaC_Accelerator::contributeInstructions(const std::string& custom_json_config) 
    {
        if (!m_isPulse || custom_json_config.empty())
        {
            return;
        }

        auto provider = xacc::getIRProvider("quantum");
        auto j = nlohmann::json::parse(custom_json_config);
        auto backends = j["backends"];
        for (auto it = backends.begin(); it != backends.end(); ++it) 
        {
            // Get the pulse library
            auto pulse_library = (*it)["specificConfiguration"]["defaults"]["pulse_library"];
            int counter = 0;
            for (auto pulse_iter = pulse_library.begin(); pulse_iter != pulse_library.end(); ++pulse_iter) 
            {
                auto pulse_name = (*pulse_iter)["name"].get<std::string>();
                auto samples = (*pulse_iter)["samples"].get<std::vector<std::vector<double>>>();

                // std::cout << counter << ", Pulse: " << pulse_name << " number of samples = " << samples.size() << " \n";

                auto pulse = std::make_shared<xacc::quantum::Pulse>(pulse_name);
                pulse->setSamples(samples);
                xacc::contributeService(pulse_name, pulse);
                counter++;
                
                {
                    const auto pulseSamples = QuaC::PulseVisitor::PulseSamplesToComplexVec(samples);
                    if (!pulseSamples.empty())
                    {
                        const auto result = m_importedPulses.emplace(pulse_name, std::move(pulseSamples));
                        if (!result.second)
                        {
                            std::cout << "Duplicate pulse with the same name '" << pulse_name << "'. The new one will be ignored.\n";
                        }
                    }
                }                
            }
        }
    }
}

