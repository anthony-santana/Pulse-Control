#pragma once

// Wrapper of the pulse provider (C++ impl)
typedef struct XaccPulseChannelProvider PulseChannelProvider;

// Get the Pulse value for a channel (unique index) at a specific time
double GetPulseValue(PulseChannelProvider* provider, int channelIdx, double time);
