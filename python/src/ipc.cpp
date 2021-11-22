#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>

#include <ipc/ipc.hpp>

namespace py = pybind11;
using namespace ipc;

bool default_can_collide(size_t, size_t) { return true; }

void define_ipc_functions(py::module_& m)
{
    m.def(
        "construct_constraint_set",
        [](const Eigen::MatrixXd& V_rest, const Eigen::MatrixXd& V,
           const Eigen::MatrixXi& E, const Eigen::MatrixXi& F, double dhat,
           const Eigen::MatrixXi& F2E, double dmin,
           const BroadPhaseMethod& method, bool ignore_codimensional_vertices,
           const std::function<bool(size_t, size_t)>& can_collide) {
            Constraints constraint_set;
            construct_constraint_set(
                V_rest, V, E, F, dhat, constraint_set, F2E, dmin, method,
                ignore_codimensional_vertices, can_collide);
            return constraint_set;
        },
        R"ipc_Qu8mg5v7(
        Construct a set of constraints used to compute the barrier potential.

        Parameters
        ----------
        V : Vertex positions as rows of a matrix
        E : Edges as rows of indicies into V
        F : Triangular faces as rows of indicies into V
        dhat : The activation distance of the barrier
        F2E : (Optional) Map from F edges to rows of E
        dmin : (Optional) Minimum distance
        ignore_codimensional_vertices : (Optional) Ignores vertices not connected to edges.
        method : (Optional) Broad-phase method to use
        can_collide : (Optional) A function that takes two vertex IDs (row numbers in F)
                      and returns true if the vertices (and faces or edges containting the
                      edges) can collide. By default all primitives can collide with all
                      other primitives.

        Returns
        -------
        The constructed set of constraints.

        See also
        --------
        Constraints

        Notes
        -----
        The given constraint_set will be cleared.
        V can either be the surface vertices or the entire mesh vertices. The edges and face should be only for the surface elements.
        )ipc_Qu8mg5v7",
        py::arg("V_rest"), py::arg("V"), py::arg("E"), py::arg("F"),
        py::arg("dhat"), py::arg("F2E") = Eigen::MatrixXi(),
        py::arg("dmin") = 0, py::arg("method") = BroadPhaseMethod::HASH_GRID,
        py::arg("ignore_codimensional_vertices") = true,
        py::arg("can_collide") =
            std::function<bool(size_t, size_t)>(default_can_collide));

    m.def(
        "compute_barrier_potential", &compute_barrier_potential,
        R"ipc_Qu8mg5v7(
        Compute the barrier potential for a given constraint set.

        Parameters
        ----------
        V : Vertex positions as rows of a matrix
        E : Edges as rows of indicies into V
        F : Triangular faces as rows of indicies into V
        constraint_set : The set of constraints
        dhat : The activation distance of the barrier

        Returns
        -------
        The sum of all barrier potentials (not scaled by the barrier stiffness).
        )ipc_Qu8mg5v7",
        py::arg("V"), py::arg("E"), py::arg("F"), py::arg("constraint_set"),
        py::arg("dhat"));

    m.def(
        "compute_barrier_potential_gradient",
        &compute_barrier_potential_gradient,
        R"ipc_Qu8mg5v7(
        Compute the gradient of the barrier potential.

        Parameters
        ----------
        V : Vertex positions as rows of a matrix
        E : Edges as rows of indicies into V
        F : Triangular faces as rows of indicies into V
        constraint_set : The set of constraints
        dhat : The activation distance of the barrier

        Returns
        -------
        The gradient of all barrier potentials (not scaled by the barrier stiffness).
        )ipc_Qu8mg5v7",
        py::arg("V"), py::arg("E"), py::arg("F"), py::arg("constraint_set"),
        py::arg("dhat"));

    m.def(
        "compute_barrier_potential_hessian", &compute_barrier_potential_hessian,
        R"ipc_Qu8mg5v7(
        Compute the hessian of the barrier potential.

        Parameters
        ----------
        V : Vertex positions as rows of a matrix
        E : Edges as rows of indicies into V
        F : Triangular faces as rows of indicies into V
        constraint_set : The set of constraints
        dhat : The activation distance of the barrier
        project_to_psd : Make sure the hessian is positive semi-definite

        Returns
        -------
        The hessian of all barrier potentials (not scaled by the barrier stiffness).
        )ipc_Qu8mg5v7",
        py::arg("V"), py::arg("E"), py::arg("F"), py::arg("constraint_set"),
        py::arg("dhat"), py::arg("project_to_psd") = true);

    // ...

    m.def(
        "has_intersections", &has_intersections,
        R"ipc_Qu8mg5v7(

        Parameters
        ----------
        V : Vertex positions as rows of a matrix
        E : Edges as rows of indicies into V
        F : Triangular faces as rows of indicies into V
        can_collide : (Optional) A function that takes two vertex IDs (row numbers in F)
                      and returns true if the vertices (and faces or edges containting the
                      edges) can collide. By default all primitives can collide with all
                      other primitives.

        Returns
        -------
        A boolean for if the mesh has intersections.
        )ipc_Qu8mg5v7",
        py::arg("V"), py::arg("E"), py::arg("F"),
        py::arg("can_collide") =
            std::function<bool(size_t, size_t)>(default_can_collide));
}
