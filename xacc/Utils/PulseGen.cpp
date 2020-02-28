#include "PulseGen.hpp"
#include <math.h> 
#include "exprtk/exprtk.hpp"

using symbol_table_t = exprtk::symbol_table<double>;
using expression_t = exprtk::expression<double>;
using parser_t = exprtk::parser<double>;

namespace QuaC {
    std::vector<std::complex<double>> SquarePulse(size_t in_nbSamples, double in_amplitude)
    {
        const std::vector<std::complex<double>> result(in_nbSamples, in_amplitude);
        return result;
    }

    std::vector<std::complex<double>> GaussianPulse(size_t in_nbSamples, double in_sigma, double in_dt, double in_amplitude)
    {
        std::vector<std::complex<double>> result;
        result.reserve(in_nbSamples);
        for (size_t i = 0; i < in_nbSamples; ++i)
        {
            result.emplace_back(in_amplitude*std::exp(-std::pow(1.0*i*in_dt, 2) / 2.0 / std::pow(in_sigma, 2.0)));
        }
        return result;
    }

    std::vector<std::complex<double>> PulseFunc(const std::string& in_functionString, size_t in_nbSamples, double in_dt)
    {
        std::vector<std::complex<double>> result;
        result.reserve(in_nbSamples);
        expression_t expression;
        parser_t parser;
        symbol_table_t symbol_table;
        double g_time = 0.0;
        symbol_table.add_variable("t", g_time);
        expression.register_symbol_table(symbol_table);
        parser.compile(in_functionString, expression);
        
        for (size_t i = 0; i < in_nbSamples; ++i)
        {
            g_time = i * in_dt;
            result.emplace_back(expression.value());
        }
        return result;
    }
}