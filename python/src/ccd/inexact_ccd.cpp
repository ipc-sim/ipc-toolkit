#include <common.hpp>

#include <ipc/config.hpp>

#include <ipc/ccd/inexact_ccd.hpp>

namespace py = pybind11;

void define_inexact_ccd(py::module_& m)
{
#ifdef IPC_TOOLKIT_WITH_INEXACT_CCD
    using namespace ipc;

    py::class_<InexactCCD, NarrowPhaseCCD>(m, "InexactCCD")
        .def(
            py::init<const double>(),
            R"ipc_Qu8mg5v7(
            Construct a new AdditiveCCD object.

            Parameters:
                conservative_rescaling: The conservative rescaling of the time of impact.
            )ipc_Qu8mg5v7",
            py::arg("conservative_rescaling") =
                InexactCCD::DEFAULT_CONSERVATIVE_RESCALING)
        .def_readonly_static(
            "DEFAULT_CONSERVATIVE_RESCALING",
            &InexactCCD::DEFAULT_CONSERVATIVE_RESCALING,
            "The default conservative rescaling value used to avoid taking steps exactly to impact.")
        .def_readonly_static(
            "SMALL_TOI", &InexactCCD::SMALL_TOI,
            "Tolerance for small time of impact which triggers rerunning CCD without a minimum separation.")
        .def_readwrite(
            "conservative_rescaling", &InexactCCD::conservative_rescaling,
            "Conservative rescaling of the time of impact.");
#endif
}
