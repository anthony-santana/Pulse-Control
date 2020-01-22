// This is the mock configs for a 2-qubit open pulse backend
// This is adapted from https://github.com/Qiskit/qiskit-terra/blob/13bc243364553667f6410b9a2f7a315c90bb598f/qiskit/test/mock/fake_openpulse_2q.py

#pragma one
#include "PulseChannelController.hpp"

struct FakePulse2Q
{
    // Ctor   
    FakePulse2Q()
    {
        {
            backendConfig.dt = 1.3333;
            // Unit: GHz
            backendConfig.loFregs_dChannels = { 4.91996800692, 5.01996800692 };
            // U channels freqs are specified by a sum with scale factors.
            backendConfig.loFregs_uChannels = { backendConfig.loFregs_dChannels[0], backendConfig.loFregs_dChannels[1] - backendConfig.loFregs_dChannels[0] };
            // Add some dummy pulse for testing
            const std::complex<double> I(0.0, 1.0);
            backendConfig.pulseLib = {
                {"test_pulse_1", std::vector<std::complex<double>>({ 0.0, 0.1*I })},
                {"test_pulse_2", std::vector<std::complex<double>>({ 0.0, 0.1*I, I })},
                {"test_pulse_3", std::vector<std::complex<double>>({ 0.0, 0.1*I, I, 0.5 })}
            };
        }
    }

    BackendChannelConfigs backendConfig;
};