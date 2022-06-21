#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>

#include <ipc/collision_mesh.hpp>

namespace py = pybind11;
using namespace ipc;

void define_collision_mesh_class(py::module_& m)
{
    py::class_<CollisionMesh>(m, "CollisionMesh")
        .def(py::init())
        .def(py::init<
             const Eigen::MatrixXd&, const Eigen::MatrixXi&,
             const Eigen::MatrixXi&>())
        .def(py::init<
             const std::vector<bool>&, const Eigen::MatrixXd&,
             const Eigen::MatrixXi&, const Eigen::MatrixXi&>())
        .def("vertices_at_rest", &CollisionMesh::vertices_at_rest)
        .def("vertices", &CollisionMesh::vertices)
        .def("edges", &CollisionMesh::edges)
        .def("faces", &CollisionMesh::faces)
        .def("faces_to_edges", &CollisionMesh::faces_to_edges)
        .def("num_vertices", &CollisionMesh::num_vertices)
        .def("dim", &CollisionMesh::dim)
        .def("ndof", &CollisionMesh::ndof)
        .def("to_full_vertex_id", &CollisionMesh::to_full_vertex_id)
        //    .def(
        //        "to_full_dof",
        //        py::overload_cast<const Eigen::VectorXd&>(
        //            &CollisionMesh::to_full_dof))
        //    .def(
        //        "to_full_dof",
        //        py::overload_cast<const Eigen::SparseMatrix<double>&>(
        //            &CollisionMesh::to_full_dof))
        .def("point_point_adjacencies", &CollisionMesh::point_point_adjacencies)
        .def("edge_point_adjacencies", &CollisionMesh::edge_point_adjacencies)
        .def("is_point_on_boundary", &CollisionMesh::is_point_on_boundary)
        .def("point_area", &CollisionMesh::point_area)
        .def("point_areas", &CollisionMesh::point_areas)
        .def("edge_area", &CollisionMesh::edge_area)
        .def("edge_areas", &CollisionMesh::edge_areas)
        .def_static(
            "construct_is_on_surface", &CollisionMesh::construct_is_on_surface)
        .def_static(
            "build_from_full_mesh", &CollisionMesh::build_from_full_mesh)
        .def_static(
            "construct_faces_to_edges",
            &CollisionMesh::construct_faces_to_edges)
        .def_readwrite("can_collide", &CollisionMesh::can_collide);
}
