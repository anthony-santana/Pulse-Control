#include "PulseChannelController.hpp"
#include <iostream>
#include <cassert>

namespace {
    const std::complex<double> I(0.0, 1.0);
}

bool BackendChannelConfigs::hasPulseName(const std::string& in_pulseName) const
{
    return pulseLib.find(in_pulseName) != pulseLib.end();
}

void BackendChannelConfigs::addOrReplacePulse(const std::string& in_pulseName, const std::vector<std::complex<double>>& in_pulseData)
{
    pulseLib[in_pulseName] = in_pulseData;
}

double BackendChannelConfigs::getPulseDuration(const std::string& in_pulseName) const
{
    if (!hasPulseName(in_pulseName))
    {
        return 0.0;
    }

    const auto iter = pulseLib.find(in_pulseName);
    return iter->second.size() * dt;
}

PulseChannelController::PulseChannelController(const BackendChannelConfigs& in_backendConfig):
    m_isInitialized(false),
    m_context({}),
    // Just put a negative time to indicate things have not been initialized yet.
    m_currentTime(-1.0),
    // Just copy configs
    m_configs(in_backendConfig)
{
    // Allocate data arrays for all channels
    m_accumulatedFcPhases = std::vector<double>(m_configs.loFregs_dChannels.size() + m_configs.loFregs_uChannels.size(), 0.0);
    m_activePulse = std::vector<const PulseScheduleEntry*>(m_configs.loFregs_dChannels.size() + m_configs.loFregs_uChannels.size(), nullptr);
    
    // Combine the list of LO freqs (D and U channels)
    m_loFreqs = m_configs.loFregs_dChannels;
    m_loFreqs.insert(m_loFreqs.end(), m_configs.loFregs_uChannels.begin(), m_configs.loFregs_uChannels.end());
}

void PulseChannelController::Initialize(const PulseScheduleRegistry& in_pulseSchedule, const FrameChangeScheduleRegistry& in_fcSchedule)
{
    m_pulseSchedules = in_pulseSchedule;
    m_fcSchedules = in_fcSchedule;
    // Move to zero time
    Tick(0.0);
}

int PulseChannelController::GetDriveChannelId(int in_dChannelIdx) const 
{
    assert(in_dChannelIdx < m_configs.loFregs_dChannels.size());
    // We use a simple global Id scheme: all D channels then U channels 
    return in_dChannelIdx;
}

int PulseChannelController::GetControlChannelId(int in_uChannelIdx) const
{
    assert(in_uChannelIdx < m_configs.loFregs_uChannels.size());
    // U channels have their global Id after D channels
    return in_uChannelIdx + m_configs.loFregs_dChannels.size();
}

double PulseChannelController::GetLoFreq(int in_channelId) const
{
    assert(in_channelId < m_loFreqs.size());
    return m_loFreqs[in_channelId];
}

void PulseChannelController::Tick(double in_time) 
{
    if (in_time <= m_currentTime)
    {
        // Nothing to do,
        // e.g. multiple channels query their data at each time step,
        // the first one has already moved the clock forward. 
        return;
    }

    // Update time
    m_currentTime = in_time;

    // Iterate over all channels (the channel is the index within the vector)
    for (size_t i = 0; i < m_activePulse.size(); ++i)
    {
        // If there is a current one that is still active:
        const auto currentPulse = m_activePulse[i];
        
        const auto isActive = [&](const PulseScheduleEntry& in_entry) -> bool {
            // Current time is within its start and stop and the conditional is satisfied
            return (in_entry.startTime <= m_currentTime) && (in_entry.stopTime >= m_currentTime) && (in_entry.conditionalFunc(m_context));
        };

        if (currentPulse != nullptr)
        {
            if (isActive(*currentPulse))
            {
                // The pulse is still active, no need to check for new schedules
                break;
            }
        }
        
        // retrieve all pulse schedules of that channel from the configs
        auto iter = m_pulseSchedules.find(i);
        if (iter != m_pulseSchedules.end())
        {
            const auto findFirstActiveEntry = [&](const std::vector<PulseScheduleEntry>& in_scheduleEntries) -> const PulseScheduleEntry* {
                for (const auto& entry: in_scheduleEntries)
                {
                    if (isActive(entry))
                    {
                        return &entry;
                    }
                }

                // Cannot find, returns nullptr
                return nullptr;        
            };

            const PulseScheduleEntry* const foundEntry = findFirstActiveEntry(iter->second);
            
            if (foundEntry != nullptr)
            {
                // Update the active pulse for the channel
                m_activePulse[i] = foundEntry;
            }
        }
    }
    
    // Update frame-change data
    {
        const auto isActivated = [&](const PulseScheduleEntry& in_entry) -> bool {
            // Current time is within its start and stop and the conditional is satisfied
            return (in_entry.startTime <= m_currentTime) && (in_entry.stopTime >= m_currentTime) && (in_entry.conditionalFunc(m_context));
        };
        
        for (size_t i = 0; i < m_accumulatedFcPhases.size(); ++i)
        { 
            auto fcIter = m_fcSchedules.find(i);
            if (fcIter != m_fcSchedules.end())
            {
               std::vector<FrameChangeCommandEntry>& fcCommandOnChannels = fcIter->second;
               for (auto it = fcCommandOnChannels.begin(); it != fcCommandOnChannels.end(); ++it)
               {
                   if (m_currentTime >= it->startTime)
                   {
                       // The FC command is activated (could be conditional)
                       if (it->conditionalFunc(m_context))
                       {
                           // Activated, add the fc phase to the accumulated phase of this channel
                           m_accumulatedFcPhases[i] += it->phase;
                       }
                       // In any cases of the conditional, after the FC command startTime (execution time),
                       // we discard this command since it has been processed.
                       fcCommandOnChannels.erase(it);
                       break;
                   }
               }
            }
        }
    }
}


double PulseChannelController::GetPulseValue(int in_channelId, double in_time)
{
    assert(in_channelId < m_activePulse.size());
    
    if (in_time > m_currentTime)
    {
        // Everytime someone querying the pulse value, we update the internal time.
        Tick(in_time);
    }
    
       
    const auto currentPulse = m_activePulse[in_channelId];
    if (currentPulse)
    {
        const auto& pulseName = currentPulse->name;
        const auto pulseArrayIter = m_configs.pulseLib.find(pulseName);

        if (pulseArrayIter == m_configs.pulseLib.end())
        {
            std::cout << "Unknown pulse named '" << pulseName << "'\n";
            return 0.0;
        }

        const auto& pulseArray = pulseArrayIter->second;
        // Compute the array index value for sampling a pulse in pulse_array.
        const auto pulseIndexCalc = [](double currentTime, const PulseScheduleEntry& pulseSchedule, size_t pulseArrayLen) -> size_t {
            return std::min(pulseArrayLen - 1, 
                            static_cast<size_t>(std::floor(1.0 * pulseArrayLen * (currentTime - pulseSchedule.startTime)/(pulseSchedule.stopTime-pulseSchedule.startTime))));
        };

        const auto pulseDataIdx = pulseIndexCalc(in_time, *currentPulse, pulseArray.size());

        const auto rawPulseData = pulseArray[pulseDataIdx];

        // Retrieve the current accumulated FC phase of this channel. 
        const double accumulatedFcPhases = m_accumulatedFcPhases[in_channelId];
        // Determine LO freq
        const auto loFreq = GetLoFreq(in_channelId);
        // Mixing signal
        const auto outputSignal = rawPulseData * std::exp(I*accumulatedFcPhases) * exp(- I * 2.0 * M_PI * loFreq * in_time);
        // Debug
        // std::cout << "D[" << in_channelId << "](" << in_time << ") = "<< outputSignal.real() << " (raw input = " << rawPulseData << "; FC phase = " << accumulatedFcPhases << ")\n";
        return outputSignal.real();
    }

    return 0.0;
}

void PulseChannelController::UpdateMeasurement(size_t in_regId, bool in_value, double in_time)
{
    // This PulseChannelController has the concept of time,
    // hence cannot go backward in time.
    assert(in_time >= m_currentTime);
    Tick(in_time);
}