#include <common.hpp>

#include <ipc/potentials/barrier_potential.hpp>

namespace py = pybind11;
using namespace ipc;

void define_barrier_potential(py::module_& m)
{
    py::class_<BarrierPotential>(m, "BarrierPotential")
        .def(
            py::init<const double>(),
            R"ipc_Qu8mg5v7(
            Construct a barrier potential.

            Parameters:
                dhat: The activation distance of the barrier.
            )ipc_Qu8mg5v7",
            py::arg("dhat"))
        .def(
            "__call__",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&,
                const CollisionConstraints&>(
                &BarrierPotential::operator(), py::const_),
            R"ipc_Qu8mg5v7(
            Compute the barrier potential for a set of contacts.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                contacts: The set of contacts.

            Returns:
                The sum of all barrier potentials (not scaled by the barrier stiffness).
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"), py::arg("contacts"))
        .def(
            "gradient",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&,
                const CollisionConstraints&>(
                &BarrierPotential::gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the barrier potential.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                contacts: The set of contacts.

            Returns:
                The gradient of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |vertices|.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"), py::arg("contacts"))
        .def(
            "hessian",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&,
                const CollisionConstraints&, const bool>(
                &BarrierPotential::hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the barrier potential.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                contacts: The set of contacts.
                project_hessian_to_psd: Make sure the hessian is positive semi-definite.

            Returns:
                The hessian of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |vertices|x|vertices|.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"), py::arg("contacts"),
            py::arg("project_hessian_to_psd") = false)
        .def(
            "shape_derivative",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&,
                const CollisionConstraints&>(
                &BarrierPotential::shape_derivative, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the shape derivative of the potential.

            std::runtime_error If the collision constraints were not built with shape derivatives enabled.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                contacts: The set of contacts.

            Returns:
                The derivative of the force with respect to X, the rest vertices.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"), py::arg("contacts"))
        .def_property(
            "dhat", [](const BarrierPotential& self) { return self.dhat(); },
            [](BarrierPotential& self, const double dhat) {
                self.set_dhat(dhat);
            },
            "Barrier activation distance.");
}
