#include <common.hpp>

#include <ipc/friction/friction_constraints.hpp>

namespace py = pybind11;
using namespace ipc;

void define_friction_constraints(py::module_& m)
{
    py::class_<FrictionConstraints>(m, "FrictionConstraints")
        .def(py::init())
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&,
                const CollisionConstraints&, double, double, double>(
                &FrictionConstraints::build),
            py::arg("mesh"), py::arg("vertices"),
            py::arg("contact_constraints"), py::arg("dhat"),
            py::arg("barrier_stiffness"), py::arg("mu"))
        .def(
            "build",
            [](FrictionConstraints& self, const CollisionMesh& mesh,
               const Eigen::MatrixXd& vertices,
               const CollisionConstraints& contact_constraints,
               const double dhat, const double barrier_stiffness,
               const Eigen::VectorXd& mus) {
                self.build(
                    mesh, vertices, contact_constraints, dhat,
                    barrier_stiffness, mus);
            },
            "", py::arg("mesh"), py::arg("vertices"),
            py::arg("contact_constraints"), py::arg("dhat"),
            py::arg("barrier_stiffness"), py::arg("mus"))
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&,
                const CollisionConstraints&, const double, const double,
                const Eigen::VectorXd&,
                const std::function<double(double, double)>&>(
                &FrictionConstraints::build),
            py::arg("mesh"), py::arg("vertices"),
            py::arg("contact_constraints"), py::arg("dhat"),
            py::arg("barrier_stiffness"), py::arg("mus"), py::arg("blend_mu"))
        .def(
            "compute_potential", &FrictionConstraints::compute_potential,
            R"ipc_Qu8mg5v7(
            Compute the friction dissapative potential from the given velocity.

            Parameters:
                mesh: The collision mesh.
                velocity: Current vertex velocity (rowwise).
                epsv: Mollifier parameter :math:`\epsilon_v`.

            Returns:
                The friction dissapative potential.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("velocity"), py::arg("epsv"))
        .def(
            "compute_potential_gradient",
            &FrictionConstraints::compute_potential_gradient,
            R"ipc_Qu8mg5v7(
            Compute the gradient of the friction dissapative potential wrt the velocity.

            Parameters:
                mesh: The collision mesh.
                velocity: Current vertex velocity (rowwise).
                epsv: Mollifier parameter :math:`\epsilon_v`.

            Returns:
                The gradient of the friction dissapative potential wrt the velocity.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("velocity"), py::arg("epsv"))
        .def(
            "compute_potential_hessian",
            &FrictionConstraints::compute_potential_hessian,
            R"ipc_Qu8mg5v7(
            Compute the Hessian of the friction dissapative potential wrt the velocity.

            Parameters:
                mesh: The collision mesh.
                velocity: Current vertex velocity (rowwise).
                epsv: Mollifier parameter :math:`\epsilon_v`.
                project_hessian_to_psd: If true, project the Hessian to be positive semi-definite.

            Returns:
                The Hessian of the friction dissapative potential wrt the velocity.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("velocity"), py::arg("epsv"),
            py::arg("project_hessian_to_psd") = false)
        .def(
            "compute_force", &FrictionConstraints::compute_force,
            R"ipc_Qu8mg5v7(
            Compute the friction force from the given velocity.

            Parameters:
                mesh: The collision mesh.
                rest_positions: Rest positions of the vertices (rowwise)
                lagged_displacements: Previous displacements of the vertices (rowwise)
                velocities: Current displacements of the vertices (rowwise)
                dhat: Barrier activation distance.
                barrier_stiffness: Barrier stiffness.
                epsv: Mollifier parameter :math:`\epsilon_v`.
                dmin: Minimum distance to use for the barrier.
                no_mu: whether to not multiply by mu

            Returns:
                The friction force.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("rest_positions"),
            py::arg("lagged_displacements"), py::arg("velocities"),
            py::arg("dhat"), py::arg("barrier_stiffness"), py::arg("epsv"),
            py::arg("dmin") = 0, py::arg("no_mu") = false)
        .def(
            "compute_force_jacobian",
            &FrictionConstraints::compute_force_jacobian,
            R"ipc_Qu8mg5v7(
            Compute the Jacobian of the friction force wrt the velocity.

            Parameters:
                mesh: The collision mesh.
                rest_positions: Rest positions of the vertices (rowwise)
                lagged_displacements: Previous displacements of the vertices (rowwise)
                velocities: Current displacements of the vertices (rowwise)
                dhat: Barrier activation distance.
                barrier_stiffness: Barrier stiffness.
                epsv: Mollifier parameter :math:`\epsilon_v`.
                wrt: The variable to take the derivative with respect to.
                dmin: Minimum distance to use for the barrier.

            Returns:
                The Jacobian of the friction force wrt the velocity.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("rest_positions"),
            py::arg("lagged_displacements"), py::arg("velocities"),
            py::arg("dhat"), py::arg("barrier_stiffness"), py::arg("epsv"),
            py::arg("wrt"), py::arg("dmin") = 0)
        .def(
            "__len__", &FrictionConstraints::size,
            "Get the number of friction constraints.")
        .def(
            "empty", &FrictionConstraints::empty,
            "Get if the friction constraints are empty.")
        .def(
            "clear", &FrictionConstraints::clear,
            "Clear the friction constraints.")
        .def(
            "__getitem__",
            [](FrictionConstraints& self, size_t idx) -> FrictionConstraint& {
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
        .def_static(
            "default_blend_mu", &FrictionConstraints::default_blend_mu,
            py::arg("mu0"), py::arg("mu1"))
        .def_readwrite("vv_constraints", &FrictionConstraints::vv_constraints)
        .def_readwrite("ev_constraints", &FrictionConstraints::ev_constraints)
        .def_readwrite("ee_constraints", &FrictionConstraints::ee_constraints)
        .def_readwrite("fv_constraints", &FrictionConstraints::fv_constraints);
}
