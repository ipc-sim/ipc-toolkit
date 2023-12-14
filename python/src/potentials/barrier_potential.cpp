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
                &BarrierPotential::Potential::operator(), py::const_),
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
                &BarrierPotential::Potential::gradient, py::const_),
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
                &BarrierPotential::Potential::hessian, py::const_),
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
        .def(
            "__call__",
            py::overload_cast<const CollisionConstraint&, const VectorMax12d&>(
                &BarrierPotential::operator(), py::const_),
            R"ipc_Qu8mg5v7(
            Compute the potential for a single contact.

            Parameters:
                contact: The contact.
                x: The contact stencil's degrees of freedom.

            Returns:
                The potential.
            )ipc_Qu8mg5v7",
            py::arg("contact"), py::arg("x"))
        .def(
            "gradient",
            py::overload_cast<const CollisionConstraint&, const VectorMax12d&>(
                &BarrierPotential::gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the potential for a single contact.

            Parameters:
                contact: The contact.
                x: The contact stencil's degrees of freedom.

            Returns:
                The gradient of the potential.
            )ipc_Qu8mg5v7",
            py::arg("contact"), py::arg("x"))
        .def(
            "hessian",
            py::overload_cast<
                const CollisionConstraint&, const VectorMax12d&, const bool>(
                &BarrierPotential::hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the potential for a single contact.

            Parameters:
                contact: The contact.
                x: The contact stencil's degrees of freedom.

            Returns:
                The hessian of the potential.
            )ipc_Qu8mg5v7",
            py::arg("contact"), py::arg("x"),
            py::arg("project_hessian_to_psd") = false)
        .def(
            "shape_derivative",
            [](const DistanceBasedPotential& self,
               const CollisionConstraint& contact,
               const std::array<long, 4>& vertex_ids,
               const VectorMax12d& rest_positions,
               const VectorMax12d& positions) {
                std::vector<Eigen::Triplet<double>> out;
                self.shape_derivative(
                    contact, vertex_ids, rest_positions, positions, out);
                Eigen::SparseMatrix<double> out_mat(
                    rest_positions.size(), rest_positions.size());
                out_mat.setFromTriplets(out.begin(), out.end());
                return out_mat;
            },
            R"ipc_Qu8mg5v7(
            Compute the shape derivative of the potential for a single contact.

            Parameters:
                contact: The contact.
                vertex_ids: The contact stencil's vertex ids.
                rest_positions: The contact stencil's rest positions.
                positions: The contact stencil's positions.
                ,out]: out Store the triplets of the shape derivative here.
            )ipc_Qu8mg5v7",
            py::arg("contact"), py::arg("vertex_ids"),
            py::arg("rest_positions"), py::arg("positions"))
        .def_property(
            "dhat", [](const BarrierPotential& self) { return self.dhat(); },
            [](BarrierPotential& self, const double dhat) {
                self.set_dhat(dhat);
            },
            "Barrier activation distance.");
}
