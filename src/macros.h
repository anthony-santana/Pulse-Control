#pragma once 

// Define some util macros which help to declare a set of channel signal functions.
// A channel function has the following signature: double Func (double time)
// All of these free functions are assigned a channel Id integer which is used
// to query the channel signal data from a centralized controller.
#define REGISTER_DRIVE_CHANNEL(pulseChannelProvider, channelId) \
    double _DriveChannel##channelId(double time) {\
        const int m_channelId = channelId;\  
        return GetPulseValue(pulseChannelProvider, m_channelId, time);\
    }

#define REGISTER_1_DRIVE_CHANNEL(pulseChannelProvider) REGISTER_DRIVE_CHANNEL(pulseChannelProvider, 0);
#define REGISTER_2_DRIVE_CHANNELS(pulseChannelProvider) \
    REGISTER_DRIVE_CHANNEL(pulseChannelProvider, 1);\
    REGISTER_1_DRIVE_CHANNEL(pulseChannelProvider) 

#define REGISTER_3_DRIVE_CHANNELS(pulseChannelProvider) \
    REGISTER_DRIVE_CHANNEL(pulseChannelProvider, 2);\
    REGISTER_2_DRIVE_CHANNELS(pulseChannelProvider) 

#define REGISTER_4_DRIVE_CHANNELS(pulseChannelProvider) \
    REGISTER_DRIVE_CHANNEL(pulseChannelProvider, 3);\
    REGISTER_3_DRIVE_CHANNELS(pulseChannelProvider) 

#define REGISTER_5_DRIVE_CHANNELS(pulseChannelProvider) \
    REGISTER_DRIVE_CHANNEL(pulseChannelProvider, 4);\
    REGISTER_4_DRIVE_CHANNELS(pulseChannelProvider) 

#define REGISTER_6_DRIVE_CHANNELS(pulseChannelProvider) \
    REGISTER_DRIVE_CHANNEL(pulseChannelProvider, 5);\
    REGISTER_5_DRIVE_CHANNELS(pulseChannelProvider) 

#define REGISTER_7_DRIVE_CHANNELS(pulseChannelProvider) \
    REGISTER_DRIVE_CHANNEL(pulseChannelProvider, 6);\
    REGISTER_6_DRIVE_CHANNELS(pulseChannelProvider) 

#define REGISTER_8_DRIVE_CHANNELS(pulseChannelProvider) \
    REGISTER_DRIVE_CHANNEL(pulseChannelProvider, 7);\
    REGISTER_7_DRIVE_CHANNELS(pulseChannelProvider)     

// Currently, support up to 8 channels
#define REGISTER_N_DRIVE_CHANNELS(pulseChannelProvider, N)  REGISTER_##N##_DRIVE_CHANNELS(pulseChannelProvider)