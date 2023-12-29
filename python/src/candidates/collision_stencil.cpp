#include <common.hpp>

#include <ipc/candidates/collision_stencil.hpp>

namespace py = pybind11;
using namespace ipc;

void define_collision_stencil(py::module_& m)
{
    py::class_<CollisionStencil>(m, "CollisionStencil")
        .def(
            "num_vertices", &CollisionStencil::num_vertices,
            "Get the number of vertices in the collision stencil.")
        .def(
            "vertex_ids", &CollisionStencil::vertex_ids,
            R"ipc_Qu8mg5v7(
            Get the vertex IDs of the collision stencil.

            Parameters:
                edges: Collision mesh edges
                faces: Collision mesh faces

            Returns:
                The vertex IDs of the collision stencil. Size is always 4, but elements i > num_vertices() are -1.
            )ipc_Qu8mg5v7",
            py::arg("edges"), py::arg("faces"))
        .def(
            "vertices", &CollisionStencil::vertices<double>,
            R"ipc_Qu8mg5v7(
            Get the vertex attributes of the collision stencil.

            T Type of the attributes

            Parameters:
                vertices: Vertex attributes
                edges: Collision mesh edges
                faces: Collision mesh faces

            Returns:
                The vertex positions of the collision stencil. Size is always 4, but elements i > num_vertices() are NaN.
            )ipc_Qu8mg5v7",
            py::arg("vertices"), py::arg("edges"), py::arg("faces"))
        .def(
            "dof", &CollisionStencil::dof<double>,
            R"ipc_Qu8mg5v7(
            Select this stencil's DOF from the full matrix of DOF.

            T Type of the DOF

            Parameters:
                X: Full matrix of DOF (rowwise).
                edges: Collision mesh edges
                faces: Collision mesh faces

            Returns:
                This stencil's DOF.
            )ipc_Qu8mg5v7",
            py::arg("X"), py::arg("edges"), py::arg("faces"))
        .def(
            "compute_distance", &CollisionStencil::compute_distance,
            R"ipc_Qu8mg5v7(
            Compute the distance of the stencil.

            Note:
                positions can be computed as stencil.dof(vertices, edges, faces)

            Parameters:
                positions: Stencil's vertex positions.

            Returns:
                Distance of the stencil.
            )ipc_Qu8mg5v7",
            py::arg("positions"))
        .def(
            "compute_distance_gradient",
            &CollisionStencil::compute_distance_gradient,
            R"ipc_Qu8mg5v7(
            Compute the distance gradient of the stencil w.r.t. the stencil's vertex positions.

            Note:
                positions can be computed as stencil.dof(vertices, edges, faces)

            Parameters:
                positions: Stencil's vertex positions.

            Returns:
                Distance gradient of the stencil w.r.t. the stencil's vertex positions.
            )ipc_Qu8mg5v7",
            py::arg("positions"))
        .def(
            "compute_distance_hessian",
            &CollisionStencil::compute_distance_hessian,
            R"ipc_Qu8mg5v7(
            Compute the distance Hessian of the stencil w.r.t. the stencil's vertex positions.

            Note:
                positions can be computed as stencil.dof(vertices, edges, faces)

            Parameters:
                positions: Stencil's vertex positions.

            Returns:
                Distance Hessian of the stencil w.r.t. the stencil's vertex positions.
            )ipc_Qu8mg5v7",
            py::arg("positions"));
}
