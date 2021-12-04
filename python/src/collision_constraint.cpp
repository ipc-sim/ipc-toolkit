#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>

#include <ipc/collision_constraint.hpp>

namespace py = pybind11;
using namespace ipc;

void define_collision_constraint_classes(py::module_& m)
{
    py::class_<CollisionConstraint>(m, "CollisionConstraint")
        .def_readwrite(
            "minimum_distance", &CollisionConstraint::minimum_distance)
        .def(
            "vertex_indices", &CollisionConstraint::vertex_indices,
            R"ipc_Qu8mg5v7(
            Get the indices of the vertices

            Parameters:
                E: edge matrix of mesh
                F: face matrix of mesh

            Returns:
                List of vertex indices
            )ipc_Qu8mg5v7",
            py::arg("E"), py::arg("F"))
        .def(
            "compute_potential", &CollisionConstraint::compute_potential,
            py::arg("V"), py::arg("E"), py::arg("F"), py::arg("dhat"))
        .def(
            "compute_potential_gradient",
            &CollisionConstraint::compute_potential_gradient, py::arg("V"),
            py::arg("E"), py::arg("F"), py::arg("dhat"))
        .def(
            "compute_potential_hessian",
            &CollisionConstraint::compute_potential_hessian, py::arg("V"),
            py::arg("E"), py::arg("F"), py::arg("dhat"),
            py::arg("project_to_psd"));

    py::class_<
        VertexVertexConstraint, VertexVertexCandidate, CollisionConstraint>(
        m, "VertexVertexConstraint")
        .def(py::init<long, long>())
        .def_readwrite("multiplicity", &VertexVertexConstraint::multiplicity);

    py::class_<EdgeVertexConstraint, EdgeVertexCandidate, CollisionConstraint>(
        m, "EdgeVertexConstraint")
        .def(py::init<long, long>())
        .def_readwrite("multiplicity", &EdgeVertexConstraint::multiplicity);

    py::class_<EdgeEdgeConstraint, EdgeEdgeCandidate, CollisionConstraint>(
        m, "EdgeEdgeConstraint")
        .def(py::init<long, long, double>())
        .def_readwrite("eps_x", &EdgeEdgeConstraint::eps_x);

    py::class_<FaceVertexConstraint, FaceVertexCandidate, CollisionConstraint>(
        m, "FaceVertexConstraint")
        .def(py::init<long, long>());

    py::class_<Constraints>(m, "Constraints")
        .def(py::init())
        .def("size", &Constraints::size)
        .def("num_constraints", &Constraints::num_constraints)
        .def("empty", &Constraints::empty)
        .def("clear", &Constraints::clear)
        .def(
            "__getitem__",
            [](Constraints& self, size_t idx) -> CollisionConstraint* {
                return &self[idx];
            })
        .def_readwrite("vv_constraints", &Constraints::vv_constraints)
        .def_readwrite("ev_constraints", &Constraints::ev_constraints)
        .def_readwrite("ee_constraints", &Constraints::ee_constraints)
        .def_readwrite("fv_constraints", &Constraints::fv_constraints);
}
