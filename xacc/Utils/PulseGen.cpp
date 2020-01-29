#include "PulseGen.hpp"

namespace QuaC {
    std::vector<std::complex<double>> SquarePulse(size_t in_nbSamples, double in_amplitude)
    {
        const std::vector<std::complex<double>> result(in_nbSamples, in_amplitude);
        return result;
    }
}