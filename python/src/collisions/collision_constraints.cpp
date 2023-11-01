#include <common.hpp>

#include <ipc/collisions/collision_constraints.hpp>

namespace py = pybind11;
using namespace ipc;

void define_collision_constraints(py::module_& m)
{
    py::class_<CollisionConstraints>(m, "CollisionConstraints")
        .def(py::init())
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&, const double,
                const double, const BroadPhaseMethod>(
                &CollisionConstraints::build),
            R"ipc_Qu8mg5v7(
            Initialize the set of constraints used to compute the barrier potential.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                dhat: The activation distance of the barrier.
                dmin: Minimum distance.
                broad_phase_method: Broad-phase method to use.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"), py::arg("dhat"),
            py::arg("dmin") = 0,
            py::arg("broad_phase_method") = DEFAULT_BROAD_PHASE_METHOD)
        .def(
            "build",
            py::overload_cast<
                const Candidates&, const CollisionMesh&, const Eigen::MatrixXd&,
                const double, const double>(&CollisionConstraints::build),
            R"ipc_Qu8mg5v7(
            Initialize the set of constraints used to compute the barrier potential.

            Parameters:
                candidates: Distance candidates from which the constraint set is built.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                dhat: The activation distance of the barrier.
                dmin:  Minimum distance.
            )ipc_Qu8mg5v7",
            py::arg("candidates"), py::arg("mesh"), py::arg("vertices"),
            py::arg("dhat"), py::arg("dmin") = 0)
        .def(
            "compute_potential", &CollisionConstraints::compute_potential,
            R"ipc_Qu8mg5v7(
            Compute the barrier potential for a given constraint set.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                dhat: The activation distance of the barrier.

            Returns:
                The sum of all barrier potentials (not scaled by the barrier stiffness).
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"), py::arg("dhat"))
        .def(
            "compute_potential_gradient",
            &CollisionConstraints::compute_potential_gradient,
            R"ipc_Qu8mg5v7(
            Compute the gradient of the barrier potential.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                dhat: The activation distance of the barrier.

            Returns:
                The gradient of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |vertices|.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"), py::arg("dhat"))
        .def(
            "compute_potential_hessian",
            &CollisionConstraints::compute_potential_hessian,
            R"ipc_Qu8mg5v7(
            Compute the hessian of the barrier potential.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                dhat: The activation distance of the barrier.
                project_hessian_to_psd: Make sure the hessian is positive semi-definite.

            Returns:
                The hessian of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |vertices|x|vertices|.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"), py::arg("dhat"),
            py::arg("project_hessian_to_psd") = false)
        .def(
            "compute_shape_derivative",
            &CollisionConstraints::compute_shape_derivative,
            R"ipc_Qu8mg5v7(
            Compute the barrier shape derivative.

            std::runtime_error If the collision constraints were not built with shape derivatives enabled.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                dhat: The activation distance of the barrier.

            Returns:
                The derivative of the force with respect to X, the rest vertices.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"), py::arg("dhat"))
        .def(
            "compute_minimum_distance",
            &CollisionConstraints::compute_minimum_distance,
            R"ipc_Qu8mg5v7(
            Computes the minimum distance between any non-adjacent elements.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The minimum distance between any non-adjacent elements.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"))
        .def(
            "__len__", &CollisionConstraints::size,
            "Get the number of collision constraints.")
        .def(
            "empty", &CollisionConstraints::empty,
            "Get if the collision constraints are empty.")
        .def(
            "clear", &CollisionConstraints::clear,
            "Clear the collision constraints.")
        .def(
            "__getitem__",
            [](CollisionConstraints& self, size_t idx) -> CollisionConstraint& {
                return self[idx];
            },
            py::return_value_policy::reference,
            R"ipc_Qu8mg5v7(
            Get a reference to constriant idx.

            Parameters:
                idx: The index of the constraint.

            Returns:
                A reference to the constraint.
            )ipc_Qu8mg5v7",
            py::arg("idx"))
        .def(
            "is_vertex_vertex", &CollisionConstraints::is_vertex_vertex,
            R"ipc_Qu8mg5v7(
            Get if the constraint at idx is a vertex-vertex constraint.

            Parameters:
                idx: The index of the constraint.

            Returns:
                If the constraint at idx is a vertex-vertex constraint.
            )ipc_Qu8mg5v7",
            py::arg("idx"))
        .def(
            "is_edge_vertex", &CollisionConstraints::is_edge_vertex,
            R"ipc_Qu8mg5v7(
            Get if the constraint at idx is an edge-vertex constraint.

            Parameters:
                idx: The index of the constraint.

            Returns:
                If the constraint at idx is an edge-vertex constraint.
            )ipc_Qu8mg5v7",
            py::arg("idx"))
        .def(
            "is_edge_edge", &CollisionConstraints::is_edge_edge,
            R"ipc_Qu8mg5v7(
            Get if the constraint at idx is an edge-edge constraint.

            Parameters:
                idx: The index of the constraint.

            Returns:
                If the constraint at idx is an edge-edge constraint.
            )ipc_Qu8mg5v7",
            py::arg("idx"))
        .def(
            "is_face_vertex", &CollisionConstraints::is_face_vertex,
            R"ipc_Qu8mg5v7(
            Get if the constraint at idx is an face-vertex constraint.

            Parameters:
                idx: The index of the constraint.

            Returns:
                If the constraint at idx is an face-vertex constraint.
            )ipc_Qu8mg5v7",
            py::arg("idx"))
        .def(
            "is_plane_vertex", &CollisionConstraints::is_plane_vertex,
            R"ipc_Qu8mg5v7(
            Get if the constraint at idx is an plane-vertex constraint.

            Parameters:
                idx: The index of the constraint.

            Returns:
                If the constraint at idx is an plane-vertex constraint.
            )ipc_Qu8mg5v7",
            py::arg("idx"))
        .def(
            "to_string", &CollisionConstraints::to_string, "", py::arg("mesh"),
            py::arg("vertices"))
        .def_property(
            "use_convergent_formulation",
            &CollisionConstraints::use_convergent_formulation,
            &CollisionConstraints::set_use_convergent_formulation,
            "If the collision constraints should use the convergent formulation.")
        .def_property(
            "are_shape_derivatives_enabled",
            &CollisionConstraints::are_shape_derivatives_enabled,
            &CollisionConstraints::set_are_shape_derivatives_enabled,
            "If the collision constraints are using the convergent formulation.")
        .def(
            "to_string", &CollisionConstraints::to_string, py::arg("mesh"),
            py::arg("vertices"))
        .def_readwrite("vv_constraints", &CollisionConstraints::vv_constraints)
        .def_readwrite("ev_constraints", &CollisionConstraints::ev_constraints)
        .def_readwrite("ee_constraints", &CollisionConstraints::ee_constraints)
        .def_readwrite("fv_constraints", &CollisionConstraints::fv_constraints)
        .def_readwrite("pv_constraints", &CollisionConstraints::pv_constraints);
}
