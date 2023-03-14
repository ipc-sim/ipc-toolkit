#include <common.hpp>

#include <ipc/candidates/candidates.hpp>

namespace py = pybind11;
using namespace ipc;

void define_candidates(py::module_& m)
{
    py::class_<Candidates>(m, "Candidates")
        .def(py::init(), "")
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&, const double,
                const BroadPhaseMethod>(&Candidates::build),
            R"ipc_Qu8mg5v7(
            Initialize the set of discrete collision detection candidates.

            Parameters:
                mesh: The surface of the contact mesh.
                vertices: Surface Vertex vertices at start as rows of a matrix.
                inflation_radius: Amount to inflate the bounding boxes.
                broad_phase_method: Broad phase method to use.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"),
            py::arg("inflation_radius") = 0,
            py::arg("broad_phase_method") = DEFAULT_BROAD_PHASE_METHOD)
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&,
                const Eigen::MatrixXd&, const double, const BroadPhaseMethod>(
                &Candidates::build),
            R"ipc_Qu8mg5v7(
            Initialize the set of continuous collision detection candidates.

            Note:
                Assumes the trajectory is linear.

            Parameters:
                mesh: The surface of the contact mesh.
                vertices_t0: Surface vertex vertices at start as rows of a matrix.
                vertices_t1: Surface vertex vertices at end as rows of a matrix.
                inflation_radius: Amount to inflate the bounding boxes.
                broad_phase_method: Broad phase method to use.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices_t0"), py::arg("vertices_t1"),
            py::arg("inflation_radius") = 0,
            py::arg("broad_phase_method") = DEFAULT_BROAD_PHASE_METHOD)
        .def("__len__", &Candidates::size, "")
        .def("empty", &Candidates::empty, "")
        .def("clear", &Candidates::clear, "")
        .def(
            "__getitem__",
            [](Candidates& self, size_t idx) -> ContinuousCollisionCandidate& {
                return self[idx];
            },
            py::return_value_policy::reference)
        .def(
            "is_step_collision_free", &Candidates::is_step_collision_free,
            R"ipc_Qu8mg5v7(
            Determine if the step is collision free from the set of candidates.

            Note:
                Assumes the trajectory is linear.

            Parameters:
                mesh: The collision mesh.
                vertices_t0: Surface vertex vertices at start as rows of a matrix.
                vertices_t1: Surface vertex vertices at end as rows of a matrix.
                min_distance: The minimum distance allowable between any two elements.
                tolerance: The tolerance for the CCD algorithm.
                max_iterations: The maximum number of iterations for the CCD algorithm.

            Returns:
                True if <b>any</b> collisions occur.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices_t0"), py::arg("vertices_t1"),
            py::arg("min_distance") = 0.0,
            py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
            py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS)
        .def(
            "compute_collision_free_stepsize",
            &Candidates::compute_collision_free_stepsize,
            R"ipc_Qu8mg5v7(
            Computes a maximal step size that is collision free using the set of collision candidates.

            Note:
                Assumes the trajectory is linear.

            Parameters:
                mesh: The collision mesh.
                vertices_t0: Vertex vertices at start as rows of a matrix. Assumes vertices_t0 is intersection free.
                vertices_t1: Surface vertex vertices at end as rows of a matrix.
                min_distance: The minimum distance allowable between any two elements.
                tolerance: The tolerance for the CCD algorithm.
                max_iterations: The maximum number of iterations for the CCD algorithm.

            Returns:
                A step-size $\in [0, 1]$ that is collision free. A value of 1.0 if a full step and 0.0 is no step.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices_t0"), py::arg("vertices_t1"),
            py::arg("min_distance") = 0.0,
            py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
            py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS)
        .def(
            "save_obj", &Candidates::save_obj, "", py::arg("filename"),
            py::arg("vertices"), py::arg("edges"), py::arg("faces"))
        .def_readwrite("ev_candidates", &Candidates::ev_candidates, "")
        .def_readwrite("ee_candidates", &Candidates::ee_candidates, "")
        .def_readwrite("fv_candidates", &Candidates::fv_candidates, "");
}
