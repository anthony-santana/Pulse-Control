#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "PulseSystemModel.hpp"
#include "xacc_service.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pyquaC, m) 
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

    py::class_<BackendChannelConfigs>(m, "BackendChannelConfigs")
        .def(py::init<>())
        .def_readwrite("loFregs_dChannels", &BackendChannelConfigs::loFregs_dChannels)
        .def_readwrite("loFregs_uChannels", &BackendChannelConfigs::loFregs_uChannels)
        .def_readwrite("dt", &BackendChannelConfigs::dt);

    py::class_<QuaC::PulseSystemModel>(m, "PulseSystemModel")
        .def(py::init<>())
        .def("name", &QuaC::PulseSystemModel::name)
        .def("setQubitT1", &QuaC::PulseSystemModel::setQubitT1)
        .def("getQubitT1", &QuaC::PulseSystemModel::getQubitT1)
        .def("loadHamiltonianJson", &QuaC::PulseSystemModel::loadHamiltonianJson)
        .def("setChannelConfigs", &QuaC::PulseSystemModel::setChannelConfigs);
}