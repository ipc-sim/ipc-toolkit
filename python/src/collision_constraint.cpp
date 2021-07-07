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
    py::class_<VertexVertexConstraint, VertexVertexCandidate>(
        m, "VertexVertexConstraint")
        .def(py::init<long, long>())
        // .def(py::init<const VertexVertexCandidate&>())
        .def(
            "vertex_indices", &VertexVertexConstraint::vertex_indices,
            R"ipc_Qu8mg5v7(
            Get the indices of the vertices

            Parameters
            ----------
            E : Edge matrix of mesh
            F : Face matrix of mesh

            Returns
            -------
            List of vertex indices
            )ipc_Qu8mg5v7",
            py::arg("E"), py::arg("F"))
        .def(
            "compute_potential", &VertexVertexConstraint::compute_potential,
            py::arg("V"), py::arg("E"), py::arg("F"), py::arg("dhat"))
        .def(
            "compute_potential_gradient",
            &VertexVertexConstraint::compute_potential_gradient, py::arg("V"),
            py::arg("E"), py::arg("F"), py::arg("dhat"))
        .def(
            "compute_potential_hessian",
            &VertexVertexConstraint::compute_potential_hessian, py::arg("V"),
            py::arg("E"), py::arg("F"), py::arg("dhat"),
            py::arg("project_to_psd"))
        .def_readwrite("multiplicity", &VertexVertexConstraint::multiplicity);
}
