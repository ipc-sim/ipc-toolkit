#include <common.hpp>

#include <ipc/potentials/friction_potential.hpp>

namespace py = pybind11;
using namespace ipc;

void define_friction_potential(py::module_& m)
{
    py::class_<FrictionPotential> friction_potential(m, "FrictionPotential");

    py::enum_<FrictionPotential::DiffWRT>(friction_potential, "DiffWRT")
        .value(
            "REST_POSITIONS", FrictionPotential::DiffWRT::REST_POSITIONS,
            "Differentiate w.r.t. rest positions")
        .value(
            "LAGGED_DISPLACEMENTS",
            FrictionPotential::DiffWRT::LAGGED_DISPLACEMENTS,
            "Differentiate w.r.t. lagged displacements")
        .value(
            "VELOCITIES", FrictionPotential::DiffWRT::VELOCITIES,
            "Differentiate w.r.t. current velocities")
        .export_values();

    friction_potential
        .def(
            py::init<const double>(),
            R"ipc_Qu8mg5v7(
            Construct a friction potential.

            Parameters:
                epsv: The smooth friction mollifier parameter :math:`\\epsilon_{v}`.
            )ipc_Qu8mg5v7",
            py::arg("epsv"))
        .def(
            "__call__",
            py::overload_cast<
                const FrictionCollisions&, const CollisionMesh&,
                const Eigen::MatrixXd&>(
                &FrictionPotential::Potential::operator(), py::const_),
            R"ipc_Qu8mg5v7(
            Compute the friction dissipative potential for a set of collisions.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The sum of all friction dissipative potentials.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"))
        .def(
            "gradient",
            py::overload_cast<
                const FrictionCollisions&, const CollisionMesh&,
                const Eigen::MatrixXd&>(
                &FrictionPotential::Potential::gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the friction dissipative potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The gradient of all friction dissipative potentials. This will have a size of |velocities|.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"))
        .def(
            "hessian",
            py::overload_cast<
                const FrictionCollisions&, const CollisionMesh&,
                const Eigen::MatrixXd&, const bool>(
                &FrictionPotential::Potential::hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the friction dissipative potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                project_hessian_to_psd: Make sure the hessian is positive semi-definite.

            Returns:
                The hessian of all friction dissipative potentials. This will have a size of |velocities|×|velocities|.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"),
            py::arg("project_hessian_to_psd") = false)
        .def(
            "force",
            py::overload_cast<
                const FrictionCollisions&, const CollisionMesh&,
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXd&, const BarrierPotential&, const double,
                const double, const bool>(
                &FrictionPotential::force, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the friction force for all collisions.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                rest_positions: Rest positions of the vertices (rowwise).
                lagged_displacements: Previous displacements of the vertices (rowwise).
                velocities: Current displacements of the vertices (rowwise).
                barrier_potential: Barrier potential (used for normal force magnitude).
                barrier_stiffness: Barrier stiffness (used for normal force magnitude).
                dmin: Minimum distance (used for normal force magnitude).
                no_mu: whether to not multiply by mu

            Returns:
                The friction force.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("rest_positions"),
            py::arg("lagged_displacements"), py::arg("velocities"),
            py::arg("barrier_potential"), py::arg("barrier_stiffness"),
            py::arg("dmin") = 0, py::arg("no_mu") = false)
        .def(
            "force_jacobian",
            py::overload_cast<
                const FrictionCollisions&, const CollisionMesh&,
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXd&, const BarrierPotential&, const double,
                const FrictionPotential::DiffWRT, const double>(
                &FrictionPotential::force_jacobian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the Jacobian of the friction force for all collisions.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                rest_positions: Rest positions of the vertices (rowwise).
                lagged_displacements: Previous displacements of the vertices (rowwise).
                velocities: Current displacements of the vertices (rowwise).
                barrier_potential: Barrier potential (used for normal force magnitude).
                barrier_stiffness: Barrier stiffness (used for normal force magnitude).
                wrt: The variable to take the derivative with respect to.
                dmin: Minimum distance (used for normal force magnitude).

            Returns:
                The Jacobian of the friction force wrt the velocities.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("rest_positions"),
            py::arg("lagged_displacements"), py::arg("velocities"),
            py::arg("barrier_potential"), py::arg("barrier_stiffness"),
            py::arg("wrt"), py::arg("dmin") = 0)
        .def(
            "__call__",
            py::overload_cast<const FrictionCollision&, const VectorMax12d&>(
                &FrictionPotential::operator(), py::const_),
            R"ipc_Qu8mg5v7(
            Compute the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The potential.
            )ipc_Qu8mg5v7",
            py::arg("collision"), py::arg("x"))
        .def(
            "gradient",
            py::overload_cast<const FrictionCollision&, const VectorMax12d&>(
                &FrictionPotential::gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The gradient of the potential.
            )ipc_Qu8mg5v7",
            py::arg("collision"), py::arg("x"))
        .def(
            "hessian",
            py::overload_cast<
                const FrictionCollision&, const VectorMax12d&, const bool>(
                &FrictionPotential::hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The hessian of the potential.
            )ipc_Qu8mg5v7",
            py::arg("collision"), py::arg("x"),
            py::arg("project_hessian_to_psd") = false)
        .def(
            "force",
            py::overload_cast<
                const FrictionCollision&, const VectorMax12d&,
                const VectorMax12d&, const VectorMax12d&,
                const BarrierPotential&, const double, const double,
                const bool>(&FrictionPotential::force, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the friction force for a single collision.

            Parameters:
                collision: The collision
                rest_positions: Rest positions of the vertices (rowwise).
                lagged_displacements: Previous displacements of the vertices (rowwise).
                velocities: Current displacements of the vertices (rowwise).
                barrier_potential: Barrier potential (used for normal force magnitude).
                barrier_stiffness: Barrier stiffness (used for normal force magnitude).
                dmin: Minimum distance (used for normal force magnitude).
                no_mu: Whether to not multiply by mu

            Returns:
                Friction force
            )ipc_Qu8mg5v7",
            py::arg("collision"), py::arg("rest_positions"),
            py::arg("lagged_displacements"), py::arg("velocities"),
            py::arg("barrier_potential"), py::arg("barrier_stiffness"),
            py::arg("dmin") = 0, py::arg("no_mu") = false)
        .def(
            "force_jacobian",
            py::overload_cast<
                const FrictionCollision&, const VectorMax12d&,
                const VectorMax12d&, const VectorMax12d&,
                const BarrierPotential&, const double,
                const FrictionPotential::DiffWRT, const double>(
                &FrictionPotential::force_jacobian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the friction force Jacobian.

            Parameters:
                collision: The collision
                rest_positions: Rest positions of the vertices (rowwise).
                lagged_displacements: Previous displacements of the vertices (rowwise).
                velocities: Current displacements of the vertices (rowwise).
                barrier_potential: Barrier potential (used for normal force magnitude).
                barrier_stiffness: Barrier stiffness (used for normal force magnitude).
                wrt: Variable to differentiate the friction force with respect to.
                dmin: Minimum distance (used for normal force magnitude).

            Returns:
                Friction force Jacobian
            )ipc_Qu8mg5v7",
            py::arg("collision"), py::arg("rest_positions"),
            py::arg("lagged_displacements"), py::arg("velocities"),
            py::arg("barrier_potential"), py::arg("barrier_stiffness"),
            py::arg("wrt"), py::arg("dmin") = 0)
        .def_property(
            "epsv", &FrictionPotential::epsv, &FrictionPotential::set_epsv,
            "The smooth friction mollifier parameter :math:`\\epsilon_{v}`.");
}
