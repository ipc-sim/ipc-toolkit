#include <common.hpp>

#include <ipc/potentials/potential.hpp>
#include <ipc/utils/hessian_assembler.hpp>

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
                Eigen::ConstRef<Eigen::MatrixXd>, const bool>(
                &Potential<TCollisions>::gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                X: Degrees of freedom of the collision mesh (e.g., vertices or velocities).
                in_full_dof: If true, return the gradient in full-mesh DOF (equivalent to mesh.to_full_dof(gradient), but assembled directly in full DOF when possible, avoiding the extra map).

            Returns:
                The gradient of the potential w.r.t. X. This will have a size of X.size (or mesh.full_ndof if in_full_dof).
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "X"_a, "in_full_dof"_a = false)
        .def(
            "hessian",
            py::overload_cast<
                const TCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>, const PSDProjectionMethod,
                const bool>(&Potential<TCollisions>::hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                X: Degrees of freedom of the collision mesh (e.g., vertices or velocities).
                project_hessian_to_psd: Make sure the hessian is positive semi-definite.
                in_full_dof: If true, return the Hessian in full-mesh DOF (equivalent to mesh.to_full_dof(hessian), but assembled directly in full DOF when possible, avoiding the two sparse-matrix products).

            Returns:
                The Hessian of the potential w.r.t. X. This will have a size of X.size by X.size (or mesh.full_ndof square if in_full_dof).
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "X"_a,
            "project_hessian_to_psd"_a = PSDProjectionMethod::NONE,
            "in_full_dof"_a = false)
        .def(
            "assemble_hessian", &Potential<TCollisions>::assemble_hessian,
            R"ipc_Qu8mg5v7(
            Assemble the Hessian of the potential using a custom assembler.

            Evaluates the local Hessian of every collision (in parallel) and feeds each to `assembler`, then leaves the
            result in the assembler (use its get_matrix()). Reuse one assembler across calls to also reuse its sparsity pattern.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                X: Degrees of freedom of the collision mesh (e.g., vertices or velocities).
                assembler: The assembler that accumulates the local Hessians.
                project_hessian_to_psd: Make sure the hessian is positive semi-definite.
                in_full_dof: If true, stencil vertex IDs are remapped to full-mesh vertex IDs (requires mesh.is_selection_dof_map).
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "X"_a, "assembler"_a,
            "project_hessian_to_psd"_a = PSDProjectionMethod::NONE,
            "in_full_dof"_a = false)
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
