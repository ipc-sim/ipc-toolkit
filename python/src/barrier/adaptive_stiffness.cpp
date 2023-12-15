#include <common.hpp>

#include <ipc/barrier/adaptive_stiffness.hpp>

namespace py = pybind11;
using namespace ipc;

void define_adaptive_stiffness(py::module_& m)
{
    m.def(
        "initial_barrier_stiffness",
        [](const double bbox_diagonal, const Barrier& barrier,
           const double dhat, const double average_mass,
           const Eigen::VectorXd& grad_energy,
           const Eigen::VectorXd& grad_barrier,
           const double min_barrier_stiffness_scale = 1e11,
           const double dmin = 0) {
            double max_barrier_stiffness;
            double r = initial_barrier_stiffness(
                bbox_diagonal, barrier, dhat, average_mass, grad_energy,
                grad_barrier, max_barrier_stiffness,
                min_barrier_stiffness_scale, dmin);
            return std::make_tuple(r, max_barrier_stiffness);
        },
        R"ipc_Qu8mg5v7(
        Compute an inital barrier stiffness using the barrier potential gradient.

        Parameters:
            bbox_diagonal: Length of the diagonal of the bounding box of the scene.
            barrier: Barrier function.
            dhat: Activation distance of the barrier.
            average_mass: Average mass of all bodies.
            grad_energy: Gradient of the elasticity energy function.
            grad_barrier: Gradient of the barrier potential.
            min_barrier_stiffness_scale: Scale used to premultiply the minimum barrier stiffness.
            dmin: Minimum distance between elements.

        Returns:
            Tuple of:
            The initial barrier stiffness.
            Maximum stiffness of the barrier.
        )ipc_Qu8mg5v7",
        py::arg("bbox_diagonal"), py::arg("barrier"), py::arg("dhat"),
        py::arg("average_mass"), py::arg("grad_energy"),
        py::arg("grad_barrier"), py::arg("min_barrier_stiffness_scale") = 1e11,
        py::arg("dmin") = 0);

    m.def(
        "update_barrier_stiffness", &update_barrier_stiffness,
        R"ipc_Qu8mg5v7(
        Update the barrier stiffness if the distance is decreasing and less than dhat_epsilon_scale * diag.

        Parameters:
            prev_min_distance: Previous minimum distance between elements.
            min_distance: Current minimum distance between elements.
            max_barrier_stiffness: Maximum stiffness of the barrier.
            barrier_stiffness: Current barrier stiffness.
            bbox_diagonal: Length of the diagonal of the bounding box of the scene.
            dhat_epsilon_scale: Update if distance is less than this fraction of the diagonal.
            dmin: Minimum distance between elements.

        Returns:
            The updated barrier stiffness.
        )ipc_Qu8mg5v7",
        py::arg("prev_min_distance"), py::arg("min_distance"),
        py::arg("max_barrier_stiffness"), py::arg("barrier_stiffness"),
        py::arg("bbox_diagonal"), py::arg("dhat_epsilon_scale") = 1e-9,
        py::arg("dmin") = 0);
}
