#include "PulseChannelController.hpp"

extern "C" {
  #include "PulseControllerHandle.h"
}


double GetPulseValue(PulseChannelProvider* provider, int channelIdx, double time)
{
    return reinterpret_cast<PulseChannelController*>(provider)->GetPulseValue(channelIdx, time);
}
