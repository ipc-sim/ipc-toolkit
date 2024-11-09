#include <common.hpp>

#include <ipc/barrier/barrier_force_magnitude.hpp>

namespace py = pybind11;
using namespace ipc;

void define_barrier_force_magnitude(py::module_& m)
{
    m.def(
        "barrier_force_magnitude", &barrier_force_magnitude,
        R"ipc_Qu8mg5v7(
        Compute the magnitude of the force due to a barrier.

        Parameters:
            distance_squared: The squared distance between elements.
            barrier: The barrier function.
            dhat: The activation distance of the barrier.
            barrier_stiffness: The stiffness of the barrier.
            dmin: The minimum distance offset to the barrier.

        Returns:
            The magnitude of the force.
        )ipc_Qu8mg5v7",
        py::arg("distance_squared"), py::arg("barrier"), py::arg("dhat"),
        py::arg("barrier_stiffness"), py::arg("dmin") = 0);

    m.def(
        "barrier_force_magnitude_gradient", &barrier_force_magnitude_gradient,
        R"ipc_Qu8mg5v7(
        Compute the gradient of the magnitude of the force due to a barrier.

        Parameters:
            distance_squared: The squared distance between elements.
            distance_squared_gradient: The gradient of the squared distance.
            barrier: The barrier function.
            dhat: The activation distance of the barrier.
            barrier_stiffness: The stiffness of the barrier.
            dmin: The minimum distance offset to the barrier.

        Returns:
            The gradient of the force.
        )ipc_Qu8mg5v7",
        py::arg("distance_squared"), py::arg("distance_squared_gradient"),
        py::arg("barrier"), py::arg("dhat"), py::arg("barrier_stiffness"),
        py::arg("dmin") = 0);
}
