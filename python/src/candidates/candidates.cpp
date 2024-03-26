#include <common.hpp>

#include <ipc/candidates/candidates.hpp>

namespace py = pybind11;
using namespace ipc;

void define_candidates(py::module_& m)
{
    py::class_<Candidates>(m, "Candidates")
        .def(py::init())
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&, const double,
                const BroadPhaseMethod>(&Candidates::build),
            R"ipc_Qu8mg5v7(
            Initialize the set of discrete collision detection candidates.

            Parameters:
                mesh: The surface of the collision mesh.
                vertices: Surface vertex positions (rowwise).
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
                mesh: The surface of the collision mesh.
                vertices_t0: Surface vertex starting positions (rowwise).
                vertices_t1: Surface vertex ending positions (rowwise).
                inflation_radius: Amount to inflate the bounding boxes.
                broad_phase_method: Broad phase method to use.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices_t0"), py::arg("vertices_t1"),
            py::arg("inflation_radius") = 0,
            py::arg("broad_phase_method") = DEFAULT_BROAD_PHASE_METHOD)
        .def("__len__", &Candidates::size)
        .def("empty", &Candidates::empty)
        .def("clear", &Candidates::clear)
        .def(
            "__getitem__",
            [](Candidates& self, size_t i) -> ContinuousCollisionCandidate& {
                return self[i];
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
                vertices_t0: Surface vertex starting positions (rowwise).
                vertices_t1: Surface vertex ending positions (rowwise).
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
                vertices_t0: Surface vertex starting positions (rowwise). Assumed to be intersection free.
                vertices_t1: Surface vertex ending positions (rowwise).
                min_distance: The minimum distance allowable between any two elements.
                tolerance: The tolerance for the CCD algorithm.
                max_iterations: The maximum number of iterations for the CCD algorithm.

            Returns:
                A step-size :math:`\in [0, 1]` that is collision free. A value of 1.0 if a full step and 0.0 is no step.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices_t0"), py::arg("vertices_t1"),
            py::arg("min_distance") = 0.0,
            py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
            py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS)
        .def(
            "compute_noncandidate_conservative_stepsize",
            &Candidates::compute_noncandidate_conservative_stepsize,
            R"ipc_Qu8mg5v7(
            Computes a conservative bound on the largest-feasible step size for surface primitives not in collision.

            Parameters:
                mesh: The collision mesh.
                displacements: Surface vertex displacements (rowwise).
                dhat: Barrier activation distance.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("displacements"), py::arg("dhat"))
        .def(
            "compute_cfl_stepsize", &Candidates::compute_cfl_stepsize,
            R"ipc_Qu8mg5v7(
            Computes a CFL-inspired CCD maximum step step size.

            Parameters:
                mesh: The collision mesh.
                vertices_t0: Surface vertex starting positions (rowwise).
                vertices_t1: Surface vertex ending positions (rowwise).
                dhat: Barrier activation distance.
                min_distance: The minimum distance allowable between any two elements.
                tolerance: The tolerance for the CCD algorithm.
                max_iterations: The maximum number of iterations for the CCD algorithm.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices_t0"), py::arg("vertices_t1"),
            py::arg("dhat"),
            py::arg("broad_phase_method") = DEFAULT_BROAD_PHASE_METHOD,
            py::arg("min_distance") = 0.0,
            py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
            py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS)
        .def(
            "save_obj", &Candidates::save_obj, py::arg("filename"),
            py::arg("vertices"), py::arg("edges"), py::arg("faces"))
        .def_readwrite("vv_candidates", &Candidates::vv_candidates)
        .def_readwrite("ev_candidates", &Candidates::ev_candidates)
        .def_readwrite("ee_candidates", &Candidates::ee_candidates)
        .def_readwrite("fv_candidates", &Candidates::fv_candidates);
}
