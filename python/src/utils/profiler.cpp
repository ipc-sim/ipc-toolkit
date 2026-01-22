#include <common.hpp>

#include <pybind11_json/pybind11_json.hpp>

#include <ipc/utils/profiler.hpp>

using namespace ipc;

void define_profiler(py::module_& m)
{
#ifdef IPC_TOOLKIT_WITH_PROFILER
    py::class_<Profiler>(m, "Profiler")
        .def("clear", &Profiler::clear, "Clear all profiling data")
        .def("reset", &Profiler::reset, "Reset the profiler data and scopes")
        .def("print", &Profiler::print, "Print the profiler data to the logger")
        .def_property_readonly(
            "csv",
            [](const Profiler& profiler) {
                std::ostringstream os;
                profiler.write_csv(os);
                return os.str();
            },
            "Get the profiler data in CSV format as a string")
        .def(
            "__str__",
            [](const Profiler& profiler) {
                std::ostringstream os;
                profiler.write_csv(os);
                return os.str();
            })
        .def_property_readonly(
            "data", py::overload_cast<>(&Profiler::data, py::const_),
            "Access the profiling data as a JSON object");

    m.def(
        "profiler", profiler, "Get the global profiler instance",
        py::return_value_policy::reference);
#endif
}
