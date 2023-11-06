#include <common.hpp>

#include <ipc/ipc.hpp>
#include <ipc/config.hpp>
#include <igl/edges.h>

namespace py = pybind11;
using namespace ipc;

void define_ipc(py::module_& m)
{
    m.attr("__version__") = IPC_TOOLKIT_VER;

    m.def(
        "is_step_collision_free", &is_step_collision_free,
        R"ipc_Qu8mg5v7(
        Determine if the step is collision free.

        Note:
            Assumes the trajectory is linear.

        Parameters:
            mesh: The collision mesh.
            vertices_t0: Surface vertex vertices at start as rows of a matrix.
            vertices_t1: Surface vertex vertices at end as rows of a matrix.
            broad_phase_method: The broad phase method to use.
            min_distance: The minimum distance allowable between any two elements.
            tolerance: The tolerance for the CCD algorithm.
            max_iterations: The maximum number of iterations for the CCD algorithm.

        Returns:
            True if <b>any</b> collisions occur.
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("vertices_t0"), py::arg("vertices_t1"),
        py::arg("broad_phase_method") = DEFAULT_BROAD_PHASE_METHOD,
        py::arg("min_distance") = 0.0,
        py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
        py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS);

    m.def(
        "compute_collision_free_stepsize", &compute_collision_free_stepsize,
        R"ipc_Qu8mg5v7(
        Computes a maximal step size that is collision free.

        Note:
            Assumes the trajectory is linear.

        Parameters:
            mesh: The collision mesh.
            vertices_t0: Vertex vertices at start as rows of a matrix. Assumes vertices_t0 is intersection free.
            vertices_t1: Surface vertex vertices at end as rows of a matrix.
            broad_phase_method: The broad phase method to use.
            min_distance: The minimum distance allowable between any two elements.
            tolerance: The tolerance for the CCD algorithm.
            max_iterations: The maximum number of iterations for the CCD algorithm.

        Returns:
            A step-size :math:`\in [0, 1]` that is collision free. A value of 1.0 if a full step and 0.0 is no step.
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("vertices_t0"), py::arg("vertices_t1"),
        py::arg("broad_phase_method") = DEFAULT_BROAD_PHASE_METHOD,
        py::arg("min_distance") = 0.0,
        py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
        py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS);

    m.def(
        "has_intersections", &has_intersections,
        R"ipc_Qu8mg5v7(
        Determine if the mesh has self intersections.

        Parameters:
            mesh: The collision mesh.
            vertices: Vertices of the collision mesh.
            broad_phase_method: The broad phase method to use.

        Returns:
            A boolean for if the mesh has intersections.
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("vertices"),
        py::arg("broad_phase_method") = DEFAULT_BROAD_PHASE_METHOD);

    m.def(
        "edges",
        [](const Eigen::MatrixXi& F) {
            Eigen::MatrixXi E;
            igl::edges(F, E);
            return E;
        },
        R"ipc_Qu8mg5v7(
        Constructs a list of unique edges represented in a given mesh F

        Parameters:
            F: #F by 3 list of mesh faces (must be triangles)

        Returns:
            #E by 2 list of edges in no particular order
        )ipc_Qu8mg5v7",
        py::arg("F"));
}
