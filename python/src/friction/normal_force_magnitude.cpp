#include <common.hpp>

#include <ipc/friction/normal_force_magnitude.hpp>

namespace py = pybind11;
using namespace ipc;

void define_normal_force_magnitude(py::module_& m)
{
    m.def(
        "compute_normal_force_magnitude", &compute_normal_force_magnitude,
        py::arg("distance_squared"), py::arg("barrier"), py::arg("dhat"),
        py::arg("barrier_stiffness"), py::arg("dmin") = 0);

    m.def(
        "compute_normal_force_magnitude_gradient",
        &compute_normal_force_magnitude_gradient, py::arg("distance_squared"),
        py::arg("distance_squared_gradient"), py::arg("barrier"),
        py::arg("dhat"), py::arg("barrier_stiffness"), py::arg("dmin") = 0);
}
