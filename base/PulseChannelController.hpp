#pragma once

#include <vector>
#include <complex>
#include <unordered_map>
#include <functional>

// Note: these data structures are designed to mimic the IBM Open Pulse specifications
// https://arxiv.org/pdf/1809.03452.pdf

// A Pulse lib is a map from string (name) to a vector of complex numbers (amplitude samples).
// The time between data points is the device time unit dt.
// Ref 5.1.5
using PulseLib = std::unordered_map<std::string, std::vector<std::complex<double>>>;

std::vector<std::complex<double>> PulseSamplesToComplexVec(const std::vector<std::vector<double>>& in_samples);

// Ref 5.1.1 (channel data only, Hamiltonian/backend device dynamic is not captured here)
struct BackendChannelConfigs
{
    // D channel LO freq.
    // Unit: GHz    
    std::vector<double> loFregs_dChannels;
    
    // U channels LO freqs 
    // Unit: GHz    
    std::vector<double> loFregs_uChannels;
    // Control signal dt (i.e. duration between pulse data points)
    // Note: this is different from the stepping dt of the solver (which should be shorter)
    // (we are simulating discrete pulse samples as provided in the QObject)  
    // Unit: ns
    double dt;   
    
    // Pulse library
    PulseLib pulseLib;

    bool hasPulseName(const std::string& in_pulseName) const;
    void addOrReplacePulse(const std::string& in_pulseName, const std::vector<std::complex<double>>& in_pulseData);
    size_t getPulseSampleSize(const std::string& in_pulseName) const;
    double getPulseDuration(const std::string& in_pulseName) const;

    template<class Archive>
    void serialize(Archive& archive)
    {
        archive(loFregs_dChannels, loFregs_uChannels, dt, pulseLib); 
    }
};


struct PulseControlContext
{
    // Register values (from measurements) which may be used to control pulse sequences.
    std::unordered_map<size_t, bool> registerValues;
    
    template<class Archive>
    void serialize(Archive& archive)
    {
        archive(registerValues); 
    }
};

// A default non-conditional entry, always returns true (i.e. activate the pulse sequence)
static const std::function<bool(const PulseControlContext&)> DEFAULT_NON_CONDITIONAL_FN = [](const PulseControlContext&){ return true; };

// A pulse scheduling entry (for a channel)
struct PulseScheduleEntry
{
    // Name of the pulse (from the library) that
    // this pulse schedule is referred to.
    std::string name;
    // Start time: when this pulse is active.
    double startTime;
    // Stop time: when to stop
    double stopTime;
    // If this function returns false, the pulse will not be activated at the specified startTime.
    std::function<bool(const PulseControlContext&)> conditionalFunc;
    // Default constructor
    PulseScheduleEntry():
     name("DEFAULT"),
     startTime(0.0),
     stopTime(0.0),
     conditionalFunc(DEFAULT_NON_CONDITIONAL_FN)
    {}

    template<class Archive>
    void serialize(Archive& archive)
    {
        archive(name, startTime, stopTime); 
    }
};

// A frame change command for a channel.
// This is a virtual command which has forward accumulated effect:
// i.e. *all* pulse data points come *after* an activated frame change command will
// have the L0 mixing phase changes.
struct FrameChangeCommandEntry
{
    // Start time: when this frame change command is active.
    double startTime;
    // Phase 
    double phase;
    // If this function returns false, the frame change will not be activated at the specified startTime.
    std::function<bool(const PulseControlContext&)> conditionalFunc;
    // Default constructor
    FrameChangeCommandEntry():
        startTime(0.0),
        phase(0.0),
        conditionalFunc(DEFAULT_NON_CONDITIONAL_FN)
    {}

    template<class Archive>
    void serialize(Archive& archive)
    {
        archive(startTime, phase); 
    }
};


using PulseScheduleRegistry = std::unordered_map<size_t, std::vector<PulseScheduleEntry>>; 
using FrameChangeScheduleRegistry = std::unordered_map<size_t, std::vector<FrameChangeCommandEntry>>; 

class PulseChannelController
{
public: 
    PulseChannelController(const BackendChannelConfigs& in_backendConfig);
    PulseChannelController(): m_currentTime(-1.0) {}
    // Initialize the pulse controller pulse/command schedules 
    void Initialize(const PulseScheduleRegistry& in_pulseSchedule, const FrameChangeScheduleRegistry& in_fcSchedule);

    // Update the pulse controller with measurement data 
    // so that it can activate/deactivate conditional pulse/command schedules appropriately.
    void UpdateMeasurement(size_t in_regId, bool in_value, double in_time);

    // This will assert that the PulseChannelController has been initialized.
    double GetPulseValue(int in_channelId, double in_time);
    // D channel
    int GetDriveChannelId(int in_dChannelIdx) const;
    // U channel
    int GetControlChannelId(int in_uChannelIdx) const;

    BackendChannelConfigs& GetBackendConfigs() { return m_configs; }

    template<class Archive>
    void serialize(Archive& archive)
    {
        archive(m_isInitialized, m_configs, m_pulseSchedules, m_fcSchedules, m_accumulatedFcPhases, m_loFreqs); 
    }

private:
    // Move the internal tracking time clock (i.e. checking all the pulse schedules w.r.t. the current time)
    void Tick(double in_time);
    double GetLoFreq(int in_channelId) const;
private:
    bool m_isInitialized;
    PulseControlContext m_context;
    double m_currentTime;
    
    // Cache of configurations
    BackendChannelConfigs m_configs;
    PulseScheduleRegistry m_pulseSchedules;
    FrameChangeScheduleRegistry m_fcSchedules;
    
    // Accumulated frame change phase of a channel
    std::vector<double> m_accumulatedFcPhases;
    // Active pulse schedule on the channel
    // This can be null if no pulse is active.
    std::vector<const PulseScheduleEntry*> m_activePulse;
    // Vector of all LO freqs (D and U channels)
    std::vector<double> m_loFreqs;
};