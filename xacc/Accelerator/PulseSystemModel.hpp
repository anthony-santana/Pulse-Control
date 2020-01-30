#pragma once

#include "Hamiltonian.hpp"
#include "PulseChannelController.hpp"
#include "CompositeInstruction.hpp"

namespace QuaC{
// A Hamiltonian model contains Hamiltonian data for pulse simulation
class HamiltonianModel
{
public:
    bool fromJson(const std::string& in_jsonString);

    std::vector<std::unique_ptr<HamiltonianTerm>>& getTerms() { return m_terms; }
private:
    std::vector<std::unique_ptr<HamiltonianTerm>> m_terms;
};


// A PulseSystemModel encapsulates all model parameters necessary for dynamical simulation.
class PulseSystemModel
{
public:
    PulseSystemModel(const std::string& in_name = "default") :
        m_name(in_name)
    {}

    // Try to load the entire system model from Json (must be in the IBM Open Pulse format)
    bool fromQobjectJson(const std::string& in_jsonString);
    
    // Try to load the Hamiltonian Json to initialize the Hamiltonian model.
    bool loadHamiltonianJson(const std::string& in_hamiltonianJsonString);
    
    // Set the channel configs 
    void setChannelConfigs(const BackendChannelConfigs& in_config);
    
    // Add a command def (as a composite instruction consisting of pulse instructions)
    // Return false if failed, e.g. the pulse instruction in the composite is not a valid one.
    bool addCommandDef(const std::string& in_cmdDefName, const std::shared_ptr<xacc::CompositeInstruction>& in_pulseComposite);

    BackendChannelConfigs& getChannelConfigs() { return m_channelConfigs; } 
    
    HamiltonianModel& getHamiltonian() { return m_hamiltonian; }
    
    // Set T1 on a qubit (this info can also be parsed from the Json)
    void setQubitT1(size_t in_qubitIdx, double in_qubitT1) { m_qubitToT1[in_qubitIdx] = in_qubitT1; }
    double getQubitT1(size_t in_qubitIdx) const;
    
private:
    std::string m_name;
    HamiltonianModel m_hamiltonian;
    BackendChannelConfigs m_channelConfigs;
    std::unordered_map<std::string, std::shared_ptr<xacc::CompositeInstruction>> m_pulseCmdDefs;
    std::unordered_map<size_t, double> m_qubitToT1;
};
}
