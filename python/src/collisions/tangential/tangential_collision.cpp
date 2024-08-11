#include <common.hpp>

#include <ipc/collisions/tangential/tangential_collision.hpp>

namespace py = pybind11;
using namespace ipc;

void define_tangential_collision(py::module_& m)
{
    py::class_<TangentialCollision, CollisionStencil>(m, "TangentialCollision")
        .def_property_readonly(
            "dim", &TangentialCollision::dim,
            "Get the dimension of the collision.")
        .def_property_readonly(
            "ndof", &TangentialCollision::ndof,
            "Get the number of degrees of freedom for the collision.")
        .def(
            "compute_tangent_basis",
            &TangentialCollision::compute_tangent_basis,
            R"ipc_Qu8mg5v7(
            Compute the tangent basis of the collision.

            Parameters:
                positions: NormalCollision stencil's vertex positions.

            Returns:
                Tangent basis of the collision.
            )ipc_Qu8mg5v7",
            py::arg("positions"))
        .def(
            "compute_tangent_basis_jacobian",
            &TangentialCollision::compute_tangent_basis_jacobian,
            R"ipc_Qu8mg5v7(
            Compute the Jacobian of the tangent basis of the collision.

            Parameters:
                positions: NormalCollision stencil's vertex positions.

            Returns:
                Jacobian of the tangent basis of the collision.
            )ipc_Qu8mg5v7",
            py::arg("positions"))
        .def(
            "compute_closest_point",
            &TangentialCollision::compute_closest_point,
            R"ipc_Qu8mg5v7(
            Compute the barycentric coordinates of the closest point.

            Parameters:
                positions: NormalCollision stencil's vertex positions.

            Returns:
                Barycentric coordinates of the closest point.
            )ipc_Qu8mg5v7",
            py::arg("positions"))
        .def(
            "compute_closest_point_jacobian",
            &TangentialCollision::compute_closest_point_jacobian,
            R"ipc_Qu8mg5v7(
            Compute the Jacobian of the barycentric coordinates of the closest point.

            Parameters:
                positions: NormalCollision stencil's vertex positions.

            Returns:
                Jacobian of the barycentric coordinates of the closest point.
            )ipc_Qu8mg5v7",
            py::arg("positions"))
        .def(
            "relative_velocity", &TangentialCollision::relative_velocity,
            R"ipc_Qu8mg5v7(
            Compute the relative velocity of the collision.

            Parameters:
                positions: NormalCollision stencil's vertex velocities.

            Returns:
                Relative velocity of the collision.
            )ipc_Qu8mg5v7",
            py::arg("velocities"))
        .def(
            "relative_velocity_matrix",
            py::overload_cast<>(
                &TangentialCollision::relative_velocity_matrix, py::const_),
            R"ipc_Qu8mg5v7(
            Construct the premultiplier matrix for the relative velocity.

            Note:
                Uses the cached closest point.

            Returns:
                A matrix M such that `relative_velocity = M * velocities`.
            )ipc_Qu8mg5v7")
        .def(
            "relative_velocity_matrix",
            py::overload_cast<const VectorMax2d&>(
                &TangentialCollision::relative_velocity_matrix, py::const_),
            R"ipc_Qu8mg5v7(
            Construct the premultiplier matrix for the relative velocity.

            Parameters:
                closest_point: Barycentric coordinates of the closest point.

            Returns:
                A matrix M such that `relative_velocity = M * velocities`.
            )ipc_Qu8mg5v7",
            py::arg("closest_point"))
        .def(
            "relative_velocity_matrix_jacobian",
            &TangentialCollision::relative_velocity_matrix_jacobian,
            R"ipc_Qu8mg5v7(
            Construct the Jacobian of the relative velocity premultiplier wrt the closest points.

            Parameters:
                closest_point: Barycentric coordinates of the closest point.

            Returns:
                Jacobian of the relative velocity premultiplier wrt the closest points.
            )ipc_Qu8mg5v7",
            py::arg("closest_point"))
        .def_readwrite(
            "normal_force_magnitude",
            &TangentialCollision::normal_force_magnitude,
            "NormalCollision force magnitude")
        .def_readwrite(
            "mu", &TangentialCollision::mu, "Coefficient of friction")
        .def_readwrite("weight", &TangentialCollision::weight, "Weight")
        .def_property(
            "weight_gradient",
            [](const TangentialCollision& self) -> Eigen::SparseMatrix<double> {
                return self.weight_gradient;
            },
            [](TangentialCollision& self,
               const Eigen::SparseMatrix<double>& weight_gradient) {
                assert_is_sparse_vector(weight_gradient, "weight_gradient");
                self.weight_gradient = weight_gradient;
            },
            "Gradient of weight with respect to all DOF")
        .def_readwrite(
            "closest_point", &TangentialCollision::closest_point,
            "Barycentric coordinates of the closest point(s)")
        .def_readwrite(
            "tangent_basis", &TangentialCollision::tangent_basis,
            "Tangent basis of the collision (max size 3×2)");
}
