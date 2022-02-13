#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>

#include <ipc/ipc.hpp>

namespace py = pybind11;
using namespace ipc;

void define_ipc_functions(py::module_& m)
{
    m.def(
        "construct_constraint_set",
        [](const CollisionMesh& mesh, const Eigen::MatrixXd& V, double dhat,
           double dmin, const BroadPhaseMethod& method) {
            Constraints constraint_set;
            construct_constraint_set(
                mesh, V, dhat, constraint_set, dmin, method);
            return constraint_set;
        },
        R"ipc_Qu8mg5v7(
        Construct a set of constraints used to compute the barrier potential.

        Parameters:
            mesh: The collision mesh.
            V: Vertices of the collision mesh.
            dhat: the activation distance of the barrier
            dmin: (optional) minimum distance
            method: (optional) broad-phase method to use

        Returns:
            The constructed set of constraints.
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("V"), py::arg("dhat"), py::arg("dmin") = 0,
        py::arg("method") = BroadPhaseMethod::HASH_GRID);

    m.def(
        "compute_barrier_potential", &compute_barrier_potential,
        R"ipc_Qu8mg5v7(
        Compute the barrier potential for a given constraint set.

        Parameters:
            mesh: The collision mesh.
            V: Vertices of the collision mesh.
            constraint_set: the set of constraints
            dhat: the activation distance of the barrier

        Returns:
            The sum of all barrier potentials (not scaled by the barrier stiffness).
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("V"), py::arg("constraint_set"),
        py::arg("dhat"));

    m.def(
        "compute_barrier_potential_gradient",
        &compute_barrier_potential_gradient,
        R"ipc_Qu8mg5v7(
        Compute the gradient of the barrier potential.

        Parameters:
            mesh: The collision mesh.
            V: Vertices of the collision mesh.
            constraint_set: the set of constraints
            dhat: the activation distance of the barrier

        Returns:
            The gradient of all barrier potentials (not scaled by the barrier stiffness).
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("V"), py::arg("constraint_set"),
        py::arg("dhat"));

    m.def(
        "compute_barrier_potential_hessian", &compute_barrier_potential_hessian,
        R"ipc_Qu8mg5v7(
        Compute the hessian of the barrier potential.

        Parameters:
            mesh: The collision mesh.
            V: Vertices of the collision mesh.
            constraint_set: the set of constraints
            dhat: the activation distance of the barrier
            project_to_psd: make sure the hessian is positive semi-definite

        Returns:
            The hessian of all barrier potentials (not scaled by the barrier stiffness).
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("V"), py::arg("constraint_set"),
        py::arg("dhat"), py::arg("project_to_psd") = true);

    // ...

    m.def(
        "has_intersections", &has_intersections,
        R"ipc_Qu8mg5v7(

        Parameters:
            mesh: The collision mesh.
            V: Vertices of the collision mesh.

        Returns:
            A boolean for if the mesh has intersections.
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("V"));
}
