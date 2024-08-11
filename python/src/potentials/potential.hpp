#include <common.hpp>

#include <ipc/potentials/potential.hpp>

namespace py = pybind11;
using namespace ipc;

template <typename TCollisions, typename PyClass>
void define_potential_methods(PyClass& potential)
{
    using TCollision = typename TCollisions::value_type;

    potential
        .def(
            "__call__",
            py::overload_cast<
                const TCollisions&, const CollisionMesh&,
                const Eigen::MatrixXd&>(
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
            py::arg("collisions"), py::arg("mesh"), py::arg("X"))
        .def(
            "gradient",
            py::overload_cast<
                const TCollisions&, const CollisionMesh&,
                const Eigen::MatrixXd&>(
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
            py::arg("collisions"), py::arg("mesh"), py::arg("X"))
        .def(
            "hessian",
            py::overload_cast<
                const TCollisions&, const CollisionMesh&,
                const Eigen::MatrixXd&, const PSDProjectionMethod>(
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
            py::arg("collisions"), py::arg("mesh"), py::arg("X"),
            py::arg("project_hessian_to_psd") = PSDProjectionMethod::NONE)
        .def(
            "__call__",
            py::overload_cast<const TCollision&, const VectorMax12d&>(
                &Potential<TCollisions>::operator(), py::const_),
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
            py::overload_cast<const TCollision&, const VectorMax12d&>(
                &Potential<TCollisions>::gradient, py::const_),
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
                const TCollision&, const VectorMax12d&,
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
            py::arg("collision"), py::arg("x"),
            py::arg("project_hessian_to_psd") = PSDProjectionMethod::NONE);
}
