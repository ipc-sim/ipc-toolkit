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
                velocities: Velocities of the collision mesh.

            Returns:
                The sum of all friction dissipative potentials.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("velocities"))
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
                velocities: Velocities of the collision mesh.

            Returns:
                The gradient of all friction dissipative potentials. This will have a size of |velocities|.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("velocities"))
        .def(
            "hessian",
            py::overload_cast<
                const FrictionCollisions&, const CollisionMesh&,
                const Eigen::MatrixXd&, const PSDProjectionMethod>(
                &FrictionPotential::Potential::hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the friction dissipative potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                velocities: Velocities of the collision mesh.
                project_hessian_to_psd: Make sure the hessian is positive semi-definite.

            Returns:
                The hessian of all friction dissipative potentials. This will have a size of |velocities|×|velocities|.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("velocities"),
            py::arg("project_hessian_to_psd") = PSDProjectionMethod::NONE)
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
                velocities: Current velocities of the vertices (rowwise).
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
}
