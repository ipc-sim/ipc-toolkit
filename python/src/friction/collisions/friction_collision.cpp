#include <common.hpp>

#include <ipc/friction/collisions/friction_collision.hpp>

namespace py = pybind11;
using namespace ipc;

void define_friction_collision(py::module_& m)
{
    py::class_<FrictionCollision, CollisionStencil>(m, "FrictionCollision")
        .def(
            "dim", &FrictionCollision::dim,
            "Get the dimension of the collision.")
        .def(
            "ndof", &FrictionCollision::ndof,
            "Get the number of degrees of freedom for the collision.")
        .def(
            "compute_normal_force_magnitude",
            &FrictionCollision::compute_normal_force_magnitude,
            R"ipc_Qu8mg5v7(
            Compute the normal force magnitude.

            Parameters:
                positions: Collision stencil's vertex positions.
                dhat: Barrier activation distance.
                barrier_stiffness: Barrier stiffness.
                dmin: Minimum distance.

            Returns:
                Normal force magnitude.
            )ipc_Qu8mg5v7",
            py::arg("positions"), py::arg("barrier_potential"),
            py::arg("barrier_stiffness"), py::arg("dmin") = 0)
        .def(
            "compute_normal_force_magnitude_gradient",
            &FrictionCollision::compute_normal_force_magnitude_gradient,
            R"ipc_Qu8mg5v7(
            Compute the gradient of the normal force magnitude.

            Parameters:
                positions: Collision stencil's vertex positions.
                dhat: Barrier activation distance.
                barrier_stiffness: Barrier stiffness.
                dmin: Minimum distance.

            Returns:
                Gradient of the normal force magnitude wrt positions.
            )ipc_Qu8mg5v7",
            py::arg("positions"), py::arg("barrier_potential"),
            py::arg("barrier_stiffness"), py::arg("dmin") = 0)
        .def(
            "compute_tangent_basis", &FrictionCollision::compute_tangent_basis,
            R"ipc_Qu8mg5v7(
            Compute the tangent basis of the collision.

            Parameters:
                positions: Collision stencil's vertex positions.

            Returns:
                Tangent basis of the collision.
            )ipc_Qu8mg5v7",
            py::arg("positions"))
        .def(
            "compute_tangent_basis_jacobian",
            &FrictionCollision::compute_tangent_basis_jacobian,
            R"ipc_Qu8mg5v7(
            Compute the Jacobian of the tangent basis of the collision.

            Parameters:
                positions: Collision stencil's vertex positions.

            Returns:
                Jacobian of the tangent basis of the collision.
            )ipc_Qu8mg5v7",
            py::arg("positions"))
        .def(
            "compute_closest_point", &FrictionCollision::compute_closest_point,
            R"ipc_Qu8mg5v7(
            Compute the barycentric coordinates of the closest point.

            Parameters:
                positions: Collision stencil's vertex positions.

            Returns:
                Barycentric coordinates of the closest point.
            )ipc_Qu8mg5v7",
            py::arg("positions"))
        .def(
            "compute_closest_point_jacobian",
            &FrictionCollision::compute_closest_point_jacobian,
            R"ipc_Qu8mg5v7(
            Compute the Jacobian of the barycentric coordinates of the closest point.

            Parameters:
                positions: Collision stencil's vertex positions.

            Returns:
                Jacobian of the barycentric coordinates of the closest point.
            )ipc_Qu8mg5v7",
            py::arg("positions"))
        .def(
            "relative_velocity", &FrictionCollision::relative_velocity,
            R"ipc_Qu8mg5v7(
            Compute the relative velocity of the collision.

            Parameters:
                positions: Collision stencil's vertex velocities.

            Returns:
                Relative velocity of the collision.
            )ipc_Qu8mg5v7",
            py::arg("velocities"))
        .def(
            "relative_velocity_matrix",
            py::overload_cast<>(
                &FrictionCollision::relative_velocity_matrix, py::const_),
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
                &FrictionCollision::relative_velocity_matrix, py::const_),
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
            &FrictionCollision::relative_velocity_matrix_jacobian,
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
            &FrictionCollision::normal_force_magnitude,
            "Collision force magnitude")
        .def_readwrite("mu", &FrictionCollision::mu, "Coefficient of friction")
        .def_readwrite("weight", &FrictionCollision::weight, "Weight")
        .def_property(
            "weight_gradient",
            [](const FrictionCollision& self) -> Eigen::SparseMatrix<double> {
                return self.weight_gradient;
            },
            [](FrictionCollision& self,
               const Eigen::SparseMatrix<double>& weight_gradient) {
                assert_is_sparse_vector(weight_gradient, "weight_gradient");
                self.weight_gradient = weight_gradient;
            },
            "Gradient of weight with respect to all DOF")
        .def_readwrite(
            "closest_point", &FrictionCollision::closest_point,
            "Barycentric coordinates of the closest point(s)")
        .def_readwrite(
            "tangent_basis", &FrictionCollision::tangent_basis,
            "Tangent basis of the collision (max size 3×2)");
}
