#include <common.hpp>

#include <ipc/ccd/additive_ccd.hpp>

namespace py = pybind11;
using namespace ipc;

void define_additive_ccd(py::module_& m)
{
    py::class_<AdditiveCCD, NarrowPhaseCCD>(m, "AdditiveCCD")
        .def(
            py::init<const double>(),
            R"ipc_Qu8mg5v7(
            Construct a new AdditiveCCD object.

            Parameters:
                conservative_rescaling: The conservative rescaling of the time of impact.
            )ipc_Qu8mg5v7",
            py::arg("conservative_rescaling") =
                AdditiveCCD::DEFAULT_CONSERVATIVE_RESCALING)
        .def_readonly_static(
            "DEFAULT_CONSERVATIVE_RESCALING",
            &AdditiveCCD::DEFAULT_CONSERVATIVE_RESCALING,
            "The default conservative rescaling value used to avoid taking steps exactly to impact. Value choosen to based on [Li et al. 2021].")
        .def_readwrite(
            "conservative_rescaling", &AdditiveCCD::conservative_rescaling,
            "The conservative rescaling value used to avoid taking steps exactly to impact.");
}
