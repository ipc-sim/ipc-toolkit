#include <common.hpp>

#include <ipc/potentials/potential.hpp>

using namespace ipc;

/// @brief Define the methods of the templated generic Potential class.
/// @tparam TCollisions Type of the collisions.
/// @tparam PyClass The pybind11 class to define the methods on.
/// @param potential The pybind11 class to define the methods on.
template <typename TCollisions, typename PyClass>
void define_potential_methods(PyClass& potential)
{
    using TCollision = typename TCollisions::value_type;

    potential
        .def(
            "__call__",
            py::overload_cast<
                const TCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>>(
                &Potential<TCollisions>::operator(), py::const_),
            R"ipc_Qu8mg5v7(
            Compute the potential for a set of collisions.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                X: Degrees of freedom of the collision mesh (e.g., vertices or velocities).

            Returns:
                The potential for a set of collisions.
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "X"_a)
        .def(
            "gradient",
            py::overload_cast<
                const TCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>>(
                &Potential<TCollisions>::gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                X: Degrees of freedom of the collision mesh (e.g., vertices or velocities).

            Returns:
                The gradient of the potential w.r.t. X. This will have a size of |X|.
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "X"_a)
        .def(
            "hessian",
            py::overload_cast<
                const TCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>, const PSDProjectionMethod>(
                &Potential<TCollisions>::hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                X: Degrees of freedom of the collision mesh (e.g., vertices or velocities).
                project_hessian_to_psd: Make sure the hessian is positive semi-definite.

            Returns:
                The Hessian of the potential w.r.t. X. This will have a size of |X|×|X|.
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "X"_a,
            "project_hessian_to_psd"_a = PSDProjectionMethod::NONE)
        .def(
            "__call__",
            py::overload_cast<const TCollision&, Eigen::ConstRef<VectorMax12d>>(
                &Potential<TCollisions>::operator(), py::const_),
            R"ipc_Qu8mg5v7(
            Compute the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The potential.
            )ipc_Qu8mg5v7",
            "collision"_a, "x"_a)
        .def(
            "gradient",
            py::overload_cast<const TCollision&, Eigen::ConstRef<VectorMax12d>>(
                &Potential<TCollisions>::gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The gradient of the potential.
            )ipc_Qu8mg5v7",
            "collision"_a, "x"_a)
        .def(
            "hessian",
            py::overload_cast<
                const TCollision&, Eigen::ConstRef<VectorMax12d>,
                const PSDProjectionMethod>(
                &Potential<TCollisions>::hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The hessian of the potential.
            )ipc_Qu8mg5v7",
            "collision"_a, "x"_a,
            "project_hessian_to_psd"_a = PSDProjectionMethod::NONE);
}
