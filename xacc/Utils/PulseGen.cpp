#include "PulseGen.hpp"
#include <math.h> 
namespace QuaC {
    std::vector<std::complex<double>> SquarePulse(size_t in_nbSamples, double in_amplitude)
    {
        const std::vector<std::complex<double>> result(in_nbSamples, in_amplitude);
        return result;
    }

    std::vector<std::complex<double>> GaussianPulse(size_t in_nbSamples, double in_sigma)
    {
        std::vector<std::complex<double>> result;
        result.reserve(in_nbSamples);
        for (size_t i = 0; i < in_nbSamples; ++i)
        {
            result.emplace_back(std::exp(-std::pow(1.0*i, 2) / 2.0 / std::pow(in_sigma, 2.0)));
        }
        return result;
    }
}