#include "common.hpp"

#include <ipc/ipc.hpp>
#include <igl/edges.h>

namespace py = pybind11;
using namespace ipc;

void define_ipc(py::module_& m)
{
    m.def(
        "compute_barrier_potential", &compute_barrier_potential,
        R"ipc_Qu8mg5v7(
        Compute the barrier potential for a given constraint set.

        Parameters:
            mesh: The collision mesh.
            V: Vertices of the collision mesh.
            constraint_set: The set of constraints.
            dhat: The activation distance of the barrier.

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
            constraint_set: The set of constraints.
            dhat: The activation distance of the barrier.

        Returns:
            The gradient of all barrier potentials (not scaled by the barrier stiffness). This will have a size of `V.size`.
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
            constraint_set: The set of constraints.
            dhat: The activation distance of the barrier.
            project_hessian_to_psd: Make sure the hessian is positive semi-definite.

        Returns:
            The hessian of all barrier potentials (not scaled by the barrier stiffness). This will have a shape of `(V.size, V.size).
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("V"), py::arg("constraint_set"),
        py::arg("dhat"), py::arg("project_hessian_to_psd") = true);

    m.def(
        "is_step_collision_free",
        py::overload_cast<
            const CollisionMesh&, const Eigen::MatrixXd&,
            const Eigen::MatrixXd&, const BroadPhaseMethod, const double,
            const long>(&is_step_collision_free),
        R"ipc_Qu8mg5v7(
        Determine if the step is collision free.

        Note:
            Assumes the trajectory is linear.

        Parameters:
            mesh: The collision mesh.
            V0: Surface vertex positions at start as rows of a matrix.
            V1: Surface vertex positions at end as rows of a matrix.

        Returns:
            True if <b>any</b> collisions occur.
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("V0"), py::arg("V1"),
        py::arg("method") = BroadPhaseMethod::HASH_GRID,
        py::arg("tolerance") = 1e-6, py::arg("max_iterations") = long(1e7));

    m.def(
        "is_step_collision_free",
        py::overload_cast<
            const Candidates&, const CollisionMesh&, const Eigen::MatrixXd&,
            const Eigen::MatrixXd&, const double, const long>(
            &is_step_collision_free),
        R"ipc_Qu8mg5v7(
        Determine if the step is collision free from a set of candidates.

        Note:
            Assumes the trajectory is linear.

        Parameters:
            candidates: Set of candidates to check for collisions.
            mesh: The collision mesh.
            V0: Surface vertex positions at start as rows of a matrix.
            V1: Surface vertex positions at end as rows of a matrix.

        Returns:
            True if <b>any</b> collisions occur.
        )ipc_Qu8mg5v7",
        py::arg("candidates"), py::arg("mesh"), py::arg("V0"), py::arg("V1"),
        py::arg("tolerance") = 1e-6, py::arg("max_iterations") = long(1e7));

    m.def(
        "compute_collision_free_stepsize",
        py::overload_cast<
            const CollisionMesh&, const Eigen::MatrixXd&,
            const Eigen::MatrixXd&, const BroadPhaseMethod, const double,
            const long>(&compute_collision_free_stepsize),
        R"ipc_Qu8mg5v7(
        Computes a maximal step size that is collision free.

        Note:
            Assumes the trajectory is linear.

        Parameters:
            mesh: The collision mesh.
            V0: Vertex positions at start as rows of a matrix. Assumes V0 is intersection free.
            V1: Surface vertex positions at end as rows of a matrix.

        Returns:
            A step-size $\in [0, 1]$ that is collision free. A value of 1.0 if a full step and 0.0 is no step.
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("V0"), py::arg("V1"),
        py::arg("method") = BroadPhaseMethod::HASH_GRID,
        py::arg("tolerance") = 1e-6, py::arg("max_iterations") = long(1e7));

    m.def(
        "compute_collision_free_stepsize",
        py::overload_cast<
            const Candidates&, const CollisionMesh&, const Eigen::MatrixXd&,
            const Eigen::MatrixXd&, const double, const long>(
            &compute_collision_free_stepsize),
        R"ipc_Qu8mg5v7(
        Computes a maximal step size that is collision free using a set of collision candidates.

        Note:
            Assumes the trajectory is linear.

        Parameters:
            candidates: Set of candidates to check for collisions.
            mesh: The collision mesh.
            V0: Vertex positions at start as rows of a matrix. Assumes V0 is intersection free.
            V1: Surface vertex positions at end as rows of a matrix.

        Returns:
            A step-size $\in [0, 1]$ that is collision free. A value of 1.0 if a full step and 0.0 is no step.
        )ipc_Qu8mg5v7",
        py::arg("candidates"), py::arg("mesh"), py::arg("V0"), py::arg("V1"),
        py::arg("tolerance") = 1e-6, py::arg("max_iterations") = long(1e7));

    m.def(
        "compute_minimum_distance", &compute_minimum_distance,
        R"ipc_Qu8mg5v7(
        Computes the minimum distance between any non-adjacent elements.

        Parameters:
            mesh: The collision mesh.
            V: Vertices of the collision mesh.

        Returns:
            The minimum distance between any non-adjacent elements.
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("V"), py::arg("constraint_set"));

    m.def(
        "has_intersections", &has_intersections,
        R"ipc_Qu8mg5v7(
        Determine if the mesh has self intersections.

        Parameters:
            mesh: The collision mesh.
            V: Vertices of the collision mesh.
            method: The broad-phase method to use.

        Returns:
            A boolean for if the mesh has intersections.
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("V"),
        py::arg("method") = BroadPhaseMethod::HASH_GRID);

    m.def(
        "edges",
        [](const Eigen::MatrixXi& F) {
            Eigen::MatrixXi E;
            igl::edges(F, E);
            return E;
        },
        R"ipc_Qu8mg5v7(
        Constructs a list of unique edges represented in a given mesh F

        Parameters:
            F: #F by 3 list of mesh faces (must be triangles)

        Returns:
            #E by 2 list of edges in no particular order
        )ipc_Qu8mg5v7",
        py::arg("F"));
}
