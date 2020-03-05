#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "PulseSystemModel.hpp"
#include "xacc_service.hpp"
#include "Pulse.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_pyquaC, m) 
{
    m.doc() = "XACC QuaC pluggin Python wrapper"; 
    m.def(
        "createPulseModel",
        []() -> std::shared_ptr<QuaC::PulseSystemModel> {
            // Get a reference to the global system model service
            auto serviceRef = xacc::getService<QuaC::PulseSystemModel>("global");
            serviceRef->reset();
            return serviceRef;
        },
    "");
    m.def(
        "addPulse",
        [](const std::string& in_pulseName, const std::vector<std::complex<double>>& in_pulseData) {
            auto pulse = std::make_shared<xacc::quantum::Pulse>(in_pulseName);
            std::vector<std::vector<double>> samples;
            samples.reserve(in_pulseData.size());
            for (const auto& dataPoint : in_pulseData)
            {
                samples.emplace_back(std::vector<double>{ dataPoint.real(), dataPoint.imag()});
            }
            pulse->setSamples(samples);
            xacc::contributeService(in_pulseName, pulse);
        },
    "");

    py::class_<BackendChannelConfigs>(m, "BackendChannelConfigs")
        .def(py::init<>())
        .def_readwrite("loFregs_dChannels", &BackendChannelConfigs::loFregs_dChannels)
        .def_readwrite("loFregs_uChannels", &BackendChannelConfigs::loFregs_uChannels)
        .def_readwrite("dt", &BackendChannelConfigs::dt);

    py::class_<QuaC::PulseSystemModel, std::shared_ptr<QuaC::PulseSystemModel>>(m, "PulseSystemModel")
        .def(py::init<>())
        .def("name", (const std::string(QuaC::PulseSystemModel::*)()) & QuaC::PulseSystemModel::name, "Get model name")
        .def("setQubitT1", (void(QuaC::PulseSystemModel::*)(size_t, double)) & QuaC::PulseSystemModel::setQubitT1, "Set T1 of a qubit")
        .def("getQubitT1", (double(QuaC::PulseSystemModel::*)(size_t)) & QuaC::PulseSystemModel::getQubitT1, "Get T1 of a qubit")
        .def("loadHamiltonianJson", (bool(QuaC::PulseSystemModel::*)(const std::string&)) & QuaC::PulseSystemModel::loadHamiltonianJson, "Load Hamiltonian from JSON")
        .def("setChannelConfigs", (void(QuaC::PulseSystemModel::*)(const BackendChannelConfigs&)) & QuaC::PulseSystemModel::setChannelConfigs, "Set backend channel configurations");
}