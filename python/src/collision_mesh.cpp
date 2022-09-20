#include "common.hpp"

#include <ipc/collision_mesh.hpp>

namespace py = pybind11;
using namespace ipc;

void define_collision_mesh(py::module_& m)
{
    py::class_<CollisionMesh>(m, "CollisionMesh")
        .def(py::init(), "")
        .def(
            py::init<
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, const Eigen::SparseMatrix<double>&>(),
            "", py::arg("vertices_at_rest"), py::arg("edges"), py::arg("faces"),
            py::arg("displacement_map") = Eigen::SparseMatrix<double>())
        .def(
            py::init<
                const std::vector<bool>&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&,
                const Eigen::SparseMatrix<double>&>(),
            "", py::arg("include_vertex"), py::arg("full_vertices_at_rest"),
            py::arg("edges"), py::arg("faces"),
            py::arg("displacement_map") = Eigen::SparseMatrix<double>())
        .def_property_readonly("num_vertices", &CollisionMesh::num_vertices, "")
        .def_property_readonly("dim", &CollisionMesh::dim, "")
        .def_property_readonly("ndof", &CollisionMesh::ndof, "")
        .def_property_readonly(
            "full_num_vertices", &CollisionMesh::full_num_vertices, "")
        .def_property_readonly("full_ndof", &CollisionMesh::full_ndof, "")
        .def_property_readonly(
            "vertices_at_rest", &CollisionMesh::vertices_at_rest, "")
        .def_property_readonly("edges", &CollisionMesh::edges, "")
        .def_property_readonly("faces", &CollisionMesh::faces, "")
        .def_property_readonly(
            "faces_to_edges", &CollisionMesh::faces_to_edges, "")
        .def("vertices", &CollisionMesh::vertices, "", py::arg("full_vertices"))
        .def(
            "displace_vertices", &CollisionMesh::displace_vertices, "",
            py::arg("full_displacements"))
        .def(
            "to_full_vertex_id", &CollisionMesh::to_full_vertex_id, "",
            py::arg("id"))
        .def(
            "to_full_dof",
            py::overload_cast<const Eigen::VectorXd&>(
                &CollisionMesh::to_full_dof, py::const_),
            "", py::arg("x"))
        .def(
            "to_full_dof",
            py::overload_cast<const Eigen::SparseMatrix<double>&>(
                &CollisionMesh::to_full_dof, py::const_),
            "", py::arg("X"))
        .def_property_readonly(
            "vertex_vertex_adjacencies",
            &CollisionMesh::vertex_vertex_adjacencies, "")
        .def_property_readonly(
            "edge_vertex_adjacencies", &CollisionMesh::edge_vertex_adjacencies,
            "")
        .def(
            "is_vertex_on_boundary", &CollisionMesh::is_vertex_on_boundary, "",
            py::arg("i"))
        .def("vertex_area", &CollisionMesh::vertex_area, "", py::arg("pi"))
        .def_property_readonly("vertex_areas", &CollisionMesh::vertex_areas, "")
        .def("edge_area", &CollisionMesh::edge_area, "", py::arg("ei"))
        .def_property_readonly("edge_areas", &CollisionMesh::edge_areas, "")
        .def_static(
            "construct_is_on_surface", &CollisionMesh::construct_is_on_surface,
            "", py::arg("num_vertices"), py::arg("edges"))
        .def_static(
            "build_from_full_mesh", &CollisionMesh::build_from_full_mesh,
            R"ipc_Qu8mg5v7(
            Helper function that automatically builds include_vertex using construct_is_on_surface.

            Parameters:
                full_vertices_at_rest: The full vertices at rest.
                edges: The edge matrix of mesh.
                faces: The face matrix of mesh.

            Returns:
                Constructed CollisionMesh.
            )ipc_Qu8mg5v7",
            py::arg("full_vertices_at_rest"), py::arg("edges"),
            py::arg("faces"))
        .def_static(
            "construct_faces_to_edges",
            &CollisionMesh::construct_faces_to_edges,
            R"ipc_Qu8mg5v7(
            Construct a matrix that maps from the faces' edges to rows in the edges matrix.

            Parameters:
                faces: The face matrix of mesh.
                edges: The edge matrix of mesh.

            Returns:
                Matrix that maps from the faces' edges to rows in the edges matrix.
            )ipc_Qu8mg5v7",
            py::arg("faces"), py::arg("edges"))
        .def_readwrite(
            "can_collide", &CollisionMesh::can_collide,
            R"ipc_Qu8mg5v7(
            A function that takes two vertex IDs (row numbers in V) and returns
            true if the vertices (and faces or edges containing the vertices) can
            collide. By default all primitives can collide with all other primitives.
            )ipc_Qu8mg5v7");
}
