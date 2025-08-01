#include <common.hpp>

#include <ipc/candidates/candidates.hpp>

using namespace ipc;

void define_candidates(py::module_& m)
{
    py::class_<Candidates>(m, "Candidates")
        .def(py::init())
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, Eigen::ConstRef<Eigen::MatrixXd>,
                const double, std::shared_ptr<BroadPhase>>(&Candidates::build),
            R"ipc_Qu8mg5v7(
            Initialize the set of discrete collision detection candidates.

            Parameters:
                mesh: The surface of the collision mesh.
                vertices: Surface vertex positions (rowwise).
                inflation_radius: Amount to inflate the bounding boxes.
                broad_phase: Broad phase to use.
            )ipc_Qu8mg5v7",
            "mesh"_a, "vertices"_a, "inflation_radius"_a = 0,
            "broad_phase"_a = make_default_broad_phase())
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXd>, const double,
                std::shared_ptr<BroadPhase>>(&Candidates::build),
            R"ipc_Qu8mg5v7(
            Initialize the set of continuous collision detection candidates.

            Note:
                Assumes the trajectory is linear.

            Parameters:
                mesh: The surface of the collision mesh.
                vertices_t0: Surface vertex starting positions (rowwise).
                vertices_t1: Surface vertex ending positions (rowwise).
                inflation_radius: Amount to inflate the bounding boxes.
                broad_phase: Broad phase to use.
            )ipc_Qu8mg5v7",
            "mesh"_a, "vertices_t0"_a, "vertices_t1"_a,
            "inflation_radius"_a = 0,
            "broad_phase"_a = make_default_broad_phase())
        .def("__len__", &Candidates::size)
        .def("empty", &Candidates::empty)
        .def("clear", &Candidates::clear)
        .def(
            "__getitem__",
            [](Candidates& self, size_t i) -> CollisionStencil& {
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
                narrow_phase_ccd: The narrow phase CCD algorithm to use.

            Returns:
                True if <b>any</b> collisions occur.
            )ipc_Qu8mg5v7",
            "mesh"_a, "vertices_t0"_a, "vertices_t1"_a, "min_distance"_a = 0.0,
            "narrow_phase_ccd"_a = DEFAULT_NARROW_PHASE_CCD)
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
                narrow_phase_ccd: The narrow phase CCD algorithm to use.

            Returns:
                A step-size :math:`\in [0, 1]` that is collision free. A value of 1.0 if a full step and 0.0 is no step.
            )ipc_Qu8mg5v7",
            "mesh"_a, "vertices_t0"_a, "vertices_t1"_a, "min_distance"_a = 0.0,
            "narrow_phase_ccd"_a = DEFAULT_NARROW_PHASE_CCD)
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
            "mesh"_a, "displacements"_a, "dhat"_a)
        .def(
            "compute_cfl_stepsize", &Candidates::compute_cfl_stepsize,
            R"ipc_Qu8mg5v7(
            Computes a CFL-inspired CCD maximum step step size.

            Parameters:
                mesh: The collision mesh.
                vertices_t0: Surface vertex starting positions (rowwise).
                vertices_t1: Surface vertex ending positions (rowwise).
                dhat: Barrier activation distance.
                min_distance: Minimum distance allowable between any two elements.
                broad_phase: Broad phase algorithm to use.
                narrow_phase_ccd: Narrow phase CCD algorithm to use.
            )ipc_Qu8mg5v7",
            "mesh"_a, "vertices_t0"_a, "vertices_t1"_a, "dhat"_a,
            "min_distance"_a = 0.0,
            "broad_phase"_a = make_default_broad_phase(),
            "narrow_phase_ccd"_a = DEFAULT_NARROW_PHASE_CCD)
        .def(
            "save_obj", &Candidates::save_obj, "filename"_a, "vertices"_a,
            "edges"_a, "faces"_a)
        .def_readwrite("vv_candidates", &Candidates::vv_candidates)
        .def_readwrite("ev_candidates", &Candidates::ev_candidates)
        .def_readwrite("ee_candidates", &Candidates::ee_candidates)
        .def_readwrite("fv_candidates", &Candidates::fv_candidates);
}
