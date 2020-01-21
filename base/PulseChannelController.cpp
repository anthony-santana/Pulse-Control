#include "PulseChannelController.hpp"
#include <iostream>
#include <cassert>

namespace {
    const std::complex<double> I(0.0, 1.0);
}

PulseChannelController::PulseChannelController():
    m_isInitialized(false),
    m_context({}),
    // Just put a negative time to indicate things have not been initialized yet.
    m_currentTime(-1.0)
{}

void PulseChannelController::Initialize(const BackendChannelConfigs& in_backendConfig, const PulseScheduleRegistry& in_pulseSchedule, const FrameChangeScheduleRegistry& in_fcSchedule)
{
    // Just copy configs
    m_configs = in_backendConfig;
    m_pulseSchedules = in_pulseSchedule;
    m_fcSchedules = in_fcSchedule;

    // TODO: we only do drive channels atm (no U channels yet)
    m_accumulatedFcPhases = std::vector<double>(m_configs.nb_dChannels, 0.0);
    m_activePulse = std::vector<const PulseScheduleEntry*>(m_configs.nb_dChannels, nullptr);
    // Move to zero time
    Tick(0.0);
}

void PulseChannelController::Tick(double in_time) 
{
    if (in_time == m_currentTime)
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
    
    // TODO: update frame-change data

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

        // TODO: we need to determine total FC phase here
        const double accumulatedFcPhases = 0.0;
        // Determine LO freq
        // TODO: need to map from channel Id to D or U channel,
        // For now, assume just D channels.
        const auto loFreq = m_configs.loFregs_dChannels[in_channelId];
        // Mixing signal
        const auto outputSignal = rawPulseData * std::exp(I*accumulatedFcPhases) * exp(- I * 2.0 * M_PI * loFreq * in_time);
        // Debug
        std::cout << "D[" << in_channelId << "](" << in_time << ") = "<< outputSignal.real() << " (raw input = " << rawPulseData << ")\n";
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