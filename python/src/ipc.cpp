#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>

#include <ipc/ipc.hpp>
#include <igl/edges.h>

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

    m.def(
        "is_step_collision_free", 
        [](const CollisionMesh& mesh, const Eigen::MatrixXd& V0,
           const Eigen::MatrixXd& V1, const BroadPhaseMethod& method,
           const double tolerance, const long max_iterations) {
        return is_step_collision_free(
            mesh, V0, V1, method, tolerance, max_iterations);
        },

        R"ipc_Qu8mg5v7(
        Determine if the step is collision free.

        V can either be the surface vertices or the entire mesh vert

        Note: Assumes the trajectory is linear.

        Parameters:
            mesh: The collision mesh.
            V0: Surface vertex positions at start as rows of a matrix.
            V1: Surface vertex positions at end as rows of a matrix.
            method: (optional) broad-phase method to use
            tolerance: (optional) query tolerance
            max_iterations: (optional) maximal number of bisection operationices.

        Returns:
            A boolean for if the mesh has intersections.
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("V0"), py::arg("V1"),
        py::arg("method") = BroadPhaseMethod::HASH_GRID,
        py::arg("tolerance") = 1e-6, py::arg("max_iterations") = 10000000);

    m.def(
        "compute_collision_free_stepsize",
        [](const CollisionMesh& mesh, const Eigen::MatrixXd& V0,
           const Eigen::MatrixXd& V1, const BroadPhaseMethod& method,
           const double tolerance, const long max_iterations) {
            return compute_collision_free_stepsize(
                mesh, V0, V1, method, tolerance, max_iterations);
        },
        R"ipc_Qu8mg5v7(
        Computes a maximal step size that is collision free.

        All vertices in V0 and V1 will be considered for collisions, so V0 and
        V1 should be only the surface vertices. The edges and face should be only
        for the surface elements.

        Note: Assumes the trajectory is linear.

        Parameters:
            mesh: The collision mesh.
            V0: Vertex positions at start as rows of a matrix. Assumes V0 is intersection free.
            V1: Surface vertex positions at end as rows of a matrix.
            method: (optional) broad-phase method to use
            tolerance: (optional) query tolerance
            max_iterations: (optional) maximal number of bisection operation

        Returns:
            A step-size \f$\in [0, 1]\f$ that is collision free. A value of 1.0 if a full step and 0.0 is no step.
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("V0"), py::arg("V1"),
        py::arg("method") = BroadPhaseMethod::HASH_GRID,
        py::arg("tolerance") = 1e-6, py::arg("max_iterations") = 10000000);

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
                py::arg("mesh"), py::arg("V"),
        py::arg("method") = BroadPhaseMethod::HASH_GRID);

    m.def(
        "edges", 
        [] (const Eigen::MatrixXi& F)
        {
            Eigen::MatrixXd E;
            igl::edges(F,E);
            return E;
        }
        ,
        R"ipc_Qu8mg5v7(
        
        Constructs a list of unique edges represented in a given mesh F

        Parameters:
            F: #F by 3 list of mesh faces (must be triangles)

        Returns:
            #E by 2 list of edges in no particular order
        )ipc_Qu8mg5v7",
                py::arg("F"));

}
        