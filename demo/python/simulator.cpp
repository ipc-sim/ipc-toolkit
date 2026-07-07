#include <ipc/dynamics/affine/joints.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>

#include <common.hpp>
#include <pybind11/stl_bind.h>
#include <pybind11_json/pybind11_json.hpp>
#include <simulator.hpp>

namespace py = pybind11;
using namespace ipc;
using namespace ipc::dynamics;
using namespace ipc::demo;

// Match the opaque declaration in dynamics/rigid/bodies.cpp so the bound
// "Poses" class is used consistently across translation units.
PYBIND11_MAKE_OPAQUE(std::vector<rigid::Pose>)

void define_simulator(py::module_& m)
{
    py::class_<Simulator> simulator(m, "Simulator", R"ipc_Qu8mg5v7(
         Body dynamics simulator with a rigid/affine toggle.

         The dynamics model is selected by Settings.body_dynamics:
         RIGID [Ferguson et al. 2021] or AFFINE [Lan et al. 2022] (which also
         supports linear equality joint constraints [Chen et al. 2022]).
         )ipc_Qu8mg5v7");

    py::enum_<Simulator::BodyDynamics>(
        simulator, "BodyDynamics",
        "The dynamics model used to simulate the bodies.")
        .value(
            "RIGID", Simulator::BodyDynamics::RIGID,
            "Rigid body dynamics [Ferguson et al. 2021]: 6-DOF per body with "
            "curved-trajectory (nonlinear) CCD.")
        .value(
            "AFFINE", Simulator::BodyDynamics::AFFINE,
            "Affine body dynamics [Lan et al. 2022]: 12-DOF per body with a "
            "stiff orthogonality potential and linear CCD.")
        .export_values();

    py::class_<Simulator::Settings>(simulator, "Settings")
        .def(py::init<>())
        .def_readwrite(
            "body_dynamics", &Simulator::Settings::body_dynamics,
            "The dynamics model used for all bodies (RIGID or AFFINE).")
        .def_readwrite(
            "bdf_order", &Simulator::Settings::bdf_order,
            "Order of the BDF time integrator (1 <= n <= 6; 1 is implicit "
            "Euler). The effective order ramps up from 1 over the first n "
            "steps.")
        .def_readwrite("dhat", &Simulator::Settings::dhat)
        .def_readwrite("gravity", &Simulator::Settings::gravity)
        .def_readwrite(
            "min_barrier_stiffness_scale",
            &Simulator::Settings::min_barrier_stiffness_scale)
        .def_readwrite(
            "dhat_epsilon_scale", &Simulator::Settings::dhat_epsilon_scale)
        .def_readwrite(
            "orthogonality_stiffness",
            &Simulator::Settings::orthogonality_stiffness,
            "Stiffness of the orthogonality potential (affine mode only).")
        .def_readwrite(
            "use_area_weighting", &Simulator::Settings::use_area_weighting)
        .def_readwrite(
            "abort_on_convergence_failure",
            &Simulator::Settings::abort_on_convergence_failure)
        .def_readwrite(
            "velocity_conv_tol", &Simulator::Settings::velocity_conv_tol,
            "Velocity convergence tolerance: converged when the proposed "
            "step moves every vertex slower than this x the bounding box "
            "diagonal (world space). Set <= 0 to disable.")
        .def_readwrite(
            "solver_params", &Simulator::Settings::solver_params,
            R"ipc_Qu8mg5v7(
             polysolve nonlinear solver parameters (a dict). Entries not set
             here get the defaults of polysolve's nonlinear-solver-spec.json.

             Note: this property returns a copy -- assign a whole dict, e.g.
             ``settings.solver_params = {**settings.solver_params, "max_iterations": 50}``.
             )ipc_Qu8mg5v7")
        .def_readwrite(
            "linear_solver_params", &Simulator::Settings::linear_solver_params,
            "polysolve linear solver parameters (a dict; empty selects "
            "polysolve's default Eigen solver). Returns a copy -- assign a "
            "whole dict.")
        .def_readwrite(
            "solver_log_level", &Simulator::Settings::solver_log_level,
            "Verbosity of polysolve's own logger (separate from the "
            "toolkit's logger; silenced by default -- lower to, e.g., "
            "ipctk.debug to see the solver's per-iteration reports).");

    simulator
        .def(
            py::init<
                const std::shared_ptr<rigid::RigidBodies>&,
                const std::vector<rigid::Pose>&, const double>(),
            "bodies"_a, "initial_poses"_a, "dt"_a)
        .def(
            py::init<
                const std::shared_ptr<rigid::RigidBodies>&,
                const std::vector<rigid::Pose>&, const double,
                const Simulator::Settings&>(),
            "bodies"_a, "initial_poses"_a, "dt"_a, "settings"_a)
        .def(
            py::init<
                const std::shared_ptr<rigid::RigidBodies>&,
                const std::vector<rigid::Pose>&,
                const std::vector<rigid::Pose>&, const double,
                const Simulator::Settings&>(),
            R"ipc_Qu8mg5v7(
             Create a simulator with initial velocities.

             Parameters:
                 bodies: Bodies in the simulation.
                 initial_poses: Initial poses of the bodies.
                 initial_velocities: Initial velocities: position is the linear velocity and rotation is the world-frame angular velocity.
                 dt: Time step.
                 settings: Simulation settings.
             )ipc_Qu8mg5v7",
            "bodies"_a, "initial_poses"_a, "initial_velocities"_a, "dt"_a,
            "settings"_a = Simulator::Settings())
        .def(
            py::init<
                const std::shared_ptr<rigid::RigidBodies>&,
                const std::vector<rigid::Pose>&,
                const std::shared_ptr<affine::JointConstraints>&, const double,
                const Simulator::Settings&>(),
            R"ipc_Qu8mg5v7(
             Create a simulator with joint constraints.

             Requires settings.body_dynamics == BodyDynamics.AFFINE
             (material-point constraints are nonlinear in the rigid rotation
             vectors).
             )ipc_Qu8mg5v7",
            "bodies"_a, "initial_poses"_a, "joints"_a, "dt"_a,
            "settings"_a = Simulator::Settings())
        .def_property("gravity", &Simulator::gravity, &Simulator::set_gravity)
        .def_property_readonly("settings", &Simulator::settings)
        .def(
            "run", &Simulator::run, "t_end"_a,
            "callback"_a = std::function<void(bool)>([](bool) { }))
        .def("step", &Simulator::step)
        .def("reset", &Simulator::reset)
        .def_property_readonly(
            "pose_history", &Simulator::pose_history,
            R"ipc_Qu8mg5v7(
             Get the history of (affine) poses in the simulation.

             In rigid mode the rotation part is exactly R(θ).

             Returns:
                 A list of affine poses at each time step.
             )ipc_Qu8mg5v7")
        .def_property_readonly(
            "rigid_pose_history", &Simulator::rigid_pose_history,
            R"ipc_Qu8mg5v7(
             Get the history of poses converted to rigid poses (log map).

             Exact in rigid mode (use for write_gltf); in affine mode the
             rotation part of the pose must be (numerically) a rotation
             matrix.

             Returns:
                 A list of poses at each time step.
             )ipc_Qu8mg5v7")
        .def_property_readonly(
            "poses", &Simulator::poses,
            R"ipc_Qu8mg5v7(
             Get the current (most recent) affine poses.

             O(num_bodies) -- prefer this over ``pose_history[-1]`` for
             per-frame access, which copies the entire history each call.
             )ipc_Qu8mg5v7")
        .def_property_readonly(
            "rigid_poses", &Simulator::rigid_poses,
            R"ipc_Qu8mg5v7(
             Get the current (most recent) poses as rigid poses (log map).

             O(num_bodies) -- the per-frame counterpart to
             ``rigid_pose_history``. Pass to ``RigidBodies.vertices`` for
             rendering.
             )ipc_Qu8mg5v7")
        .def_property_readonly("t", &Simulator::t);
}
