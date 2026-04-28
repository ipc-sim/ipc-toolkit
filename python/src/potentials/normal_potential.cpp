#include <common.hpp>
#include <potentials/potential.hpp>

#include <ipc/potentials/normal_potential.hpp>

using namespace ipc;

void define_normal_potential(py::module_& m)
{
    py::class_<NormalPotential> normal_potential(m, "NormalPotential");

    define_potential_methods<NormalCollisions>(normal_potential);

    normal_potential
        .def(
            "shape_derivative",
            py::overload_cast<
                const NormalCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>>(
                &NormalPotential::shape_derivative, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the shape derivative of the potential.

            std::runtime_error If the collision collisions were not built with shape derivatives enabled.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The derivative of the force with respect to X, the rest vertices.
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "vertices"_a)
        .def(
            "shape_derivative",
            [](const NormalPotential& self, const NormalCollision& collision,
               const std::array<index_t, 4>& vertex_ids,
               Eigen::ConstRef<VectorMax12d> rest_positions,
               Eigen::ConstRef<VectorMax12d> positions) {
                std::vector<Eigen::Triplet<double>> out;
                self.shape_derivative(
                    collision, vertex_ids, rest_positions, positions, out);
                return out;
            },
            R"ipc_Qu8mg5v7(
            Compute the shape derivative of the potential for a single collision.

            Parameters:
                collision: The collision.
                vertex_ids: The collision stencil's vertex ids.
                rest_positions: The collision stencil's rest positions.
                positions: The collision stencil's positions.
                ,out]: out Store the triplets of the shape derivative here.
            )ipc_Qu8mg5v7",
            "collision"_a, "vertex_ids"_a, "rest_positions"_a, "positions"_a)
        .def(
            "force_magnitude", &NormalPotential::force_magnitude,
            R"ipc_Qu8mg5v7(
            Compute the force magnitude for a collision.

            Parameters:
                distance_squared: The squared distance between elements.
                dmin: The minimum distance offset to the barrier.
                barrier_stiffness: The barrier stiffness.

            Returns:
                The force magnitude.
            )ipc_Qu8mg5v7",
            "distance_squared"_a, "dmin"_a)
        .def(
            "force_magnitude_gradient",
            &NormalPotential::force_magnitude_gradient,
            R"ipc_Qu8mg5v7(
            Compute the gradient of the force magnitude for a collision.

            Parameters:
                distance_squared: The squared distance between elements.
                distance_squared_gradient: The gradient of the squared distance.
                dmin: The minimum distance offset to the barrier.

            Returns:
                The gradient of the force.
            )ipc_Qu8mg5v7",
            "distance_squared"_a, "distance_squared_gradient"_a, "dmin"_a)
        .def(
            "gauss_newton_hessian_diagonal",
            py::overload_cast<
                const NormalCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>>(
                &NormalPotential::gauss_newton_hessian_diagonal, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the diagonal of the cumulative Gauss-Newton Hessian of the potential.

            Uses the distance-vector formulation to efficiently compute the
            diagonal of the Gauss-Newton Hessian without forming full local
            12x12 Hessian matrices. This is useful as a Jacobi preconditioner
            for iterative solvers.

            Note:
                This is a Gauss-Newton approximation (drops derivatives of
                closest-point coefficients), not the exact Hessian.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The diagonal of the Gauss-Newton Hessian as a vector of size vertices.size().
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "vertices"_a)
        .def(
            "gauss_newton_hessian_diagonal",
            py::overload_cast<
                const NormalCollision&, Eigen::ConstRef<VectorMax12d>>(
                &NormalPotential::gauss_newton_hessian_diagonal, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the diagonal of the Gauss-Newton Hessian for a single collision.

            Uses the distance-vector formulation (Eqs. 10-12) to efficiently
            compute the diagonal without forming the full local Hessian.

            Note:
                This is a Gauss-Newton approximation, not the exact Hessian diagonal.

            Parameters:
                collision: The collision.
                positions: The collision stencil's positions.

            Returns:
                The diagonal of the Gauss-Newton Hessian as a vector of size ndof.
            )ipc_Qu8mg5v7",
            "collision"_a, "positions"_a)
        .def(
            "gauss_newton_hessian_quadratic_form",
            py::overload_cast<
                const NormalCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::VectorXd>>(
                &NormalPotential::gauss_newton_hessian_quadratic_form,
                py::const_),
            R"ipc_Qu8mg5v7(
            Compute the product p^T H p for the cumulative Gauss-Newton Hessian.

            Uses the distance-vector formulation to efficiently compute the
            quadratic form without forming full local 12x12 Hessian matrices
            nor the global sparse Hessian. This is useful for nonlinear
            conjugate gradient methods.

            Note:
                This is a Gauss-Newton approximation (drops derivatives of
                closest-point coefficients), not the exact Hessian.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                p: The direction vector of size vertices.size().

            Returns:
                The scalar value p^T H p (approximate).
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "vertices"_a, "p"_a)
        .def(
            "gauss_newton_hessian_quadratic_form",
            py::overload_cast<
                const NormalCollision&, Eigen::ConstRef<VectorMax12d>,
                Eigen::ConstRef<VectorMax12d>>(
                &NormalPotential::gauss_newton_hessian_quadratic_form,
                py::const_),
            R"ipc_Qu8mg5v7(
            Compute p^T H p for a single collision using the Gauss-Newton Hessian.

            Uses the distance-vector formulation (Eqs. 10, 13-14) to
            efficiently compute the quadratic form without forming the full
            local Hessian.

            Note:
                This is a Gauss-Newton approximation, not the exact Hessian.

            Parameters:
                collision: The collision.
                positions: The collision stencil's positions.
                p: The local direction vector (size ndof).

            Returns:
                The scalar value p^T H p (approximate).
            )ipc_Qu8mg5v7",
            "collision"_a, "positions"_a, "p"_a);
}
