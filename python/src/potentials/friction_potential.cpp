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
            "Compute the force jacobian w.r.t. rest positions")
        .value(
            "LAGGED_DISPLACEMENTS",
            FrictionPotential::DiffWRT::LAGGED_DISPLACEMENTS,
            "Compute the force jacobian w.r.t. rest lagged displacements")
        .value(
            "VELOCITIES", FrictionPotential::DiffWRT::VELOCITIES,
            "Compute the force jacobian w.r.t. velocities")
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
            Compute the barrier potential for a set of collisions.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The sum of all barrier potentials (not scaled by the barrier stiffness).
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"))
        .def(
            "gradient",
            py::overload_cast<
                const FrictionCollisions&, const CollisionMesh&,
                const Eigen::MatrixXd&>(
                &FrictionPotential::Potential::gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the barrier potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The gradient of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |vertices|.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"))
        .def(
            "hessian",
            py::overload_cast<
                const FrictionCollisions&, const CollisionMesh&,
                const Eigen::MatrixXd&, const bool>(
                &FrictionPotential::Potential::hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the barrier potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                project_hessian_to_psd: Make sure the hessian is positive semi-definite.

            Returns:
                The hessian of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |vertices|x|vertices|.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"),
            py::arg("project_hessian_to_psd") = false)
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
        .def_property(
            "epsv", [](const FrictionPotential& self) { return self.epsv(); },
            [](FrictionPotential& self, const double epsv) {
                self.set_epsv(epsv);
            },
            "The smooth friction mollifier parameter :math:`\\epsilon_{v}`.");
}
