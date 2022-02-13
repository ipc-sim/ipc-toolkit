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
        .def("edges", &CollisionMesh::edges)
        .def("faces", &CollisionMesh::faces)
        .def("faces_to_edges", &CollisionMesh::faces_to_edges)
        .def("num_vertices", &CollisionMesh::num_vertices)
        .def("dim", &CollisionMesh::dim)
        .def("ndof", &CollisionMesh::ndof)
        .def("to_full_vertex_id", &CollisionMesh::to_full_vertex_id)
        .def_readwrite("can_collide", &CollisionMesh::can_collide)
        //
        ;
}
