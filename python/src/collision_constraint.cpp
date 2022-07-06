#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>

#include <ipc/src/collision_constraint.hpp>

namespace py = pybind11;
using namespace ipc;

void define_collision_constraint_members(py::module_& m)
{
    py::class_<CollisionConstraint>(m, "CollisionConstraint")
        .def("num_vertices", &CollisionConstraint::num_vertices, "")
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
            "compute_distance", &CollisionConstraint::compute_distance, "",
            py::arg("V"), py::arg("E"), py::arg("F"))
        .def(
            "compute_distance_gradient",
            &CollisionConstraint::compute_distance_gradient, "", py::arg("V"),
            py::arg("E"), py::arg("F"))
        .def(
            "compute_distance_hessian",
            &CollisionConstraint::compute_distance_hessian, "", py::arg("V"),
            py::arg("E"), py::arg("F"))
        .def(
            "compute_potential", &CollisionConstraint::compute_potential, "",
            py::arg("V"), py::arg("E"), py::arg("F"), py::arg("dhat"))
        .def(
            "compute_potential_gradient",
            &CollisionConstraint::compute_potential_gradient, "", py::arg("V"),
            py::arg("E"), py::arg("F"), py::arg("dhat"))
        .def(
            "compute_potential_hessian",
            &CollisionConstraint::compute_potential_hessian, "", py::arg("V"),
            py::arg("E"), py::arg("F"), py::arg("dhat"),
            py::arg("project_hessian_to_psd"))
        .def_readwrite(
            "minimum_distance", &CollisionConstraint::minimum_distance, "");

    py::class_<
        VertexVertexConstraint, VertexVertexCandidate, CollisionConstraint>(
        m, "VertexVertexConstraint")
        .def(
            py::init<long, long>(), "", py::arg("vertex0_index"),
            py::arg("vertex1_index"))
        .def(py::init<const VertexVertexCandidate&>(), "", py::arg("candidate"))
        .def("num_vertices", &VertexVertexConstraint::num_vertices, "")
        .def(
            "vertex_indices", &VertexVertexConstraint::vertex_indices, "",
            py::arg("E"), py::arg("F"))
        .def(
            "compute_distance", &VertexVertexConstraint::compute_distance, "",
            py::arg("V"), py::arg("E"), py::arg("F"))
        .def(
            "compute_distance_gradient",
            &VertexVertexConstraint::compute_distance_gradient, "",
            py::arg("V"), py::arg("E"), py::arg("F"))
        .def(
            "compute_distance_hessian",
            &VertexVertexConstraint::compute_distance_hessian, "", py::arg("V"),
            py::arg("E"), py::arg("F"))
        .def(
            "compute_potential", &VertexVertexConstraint::compute_potential, "",
            py::arg("V"), py::arg("E"), py::arg("F"), py::arg("dhat"))
        .def(
            "compute_potential_gradient",
            &VertexVertexConstraint::compute_potential_gradient, "",
            py::arg("V"), py::arg("E"), py::arg("F"), py::arg("dhat"))
        .def(
            "compute_potential_hessian",
            &VertexVertexConstraint::compute_potential_hessian, "",
            py::arg("V"), py::arg("E"), py::arg("F"), py::arg("dhat"),
            py::arg("project_hessian_to_psd"))
        .def_readwrite(
            "multiplicity", &VertexVertexConstraint::multiplicity, "");

    py::class_<EdgeVertexConstraint, EdgeVertexCandidate, CollisionConstraint>(
        m, "EdgeVertexConstraint")
        .def(
            py::init<long, long>(), "", py::arg("edge_index"),
            py::arg("vertex_index"))
        .def(py::init<const EdgeVertexCandidate&>(), "", py::arg("candidate"))
        .def("num_vertices", &EdgeVertexConstraint::num_vertices, "")
        .def(
            "vertex_indices", &EdgeVertexConstraint::vertex_indices, "",
            py::arg("E"), py::arg("F"))
        .def(
            "compute_distance", &EdgeVertexConstraint::compute_distance, "",
            py::arg("V"), py::arg("E"), py::arg("F"))
        .def(
            "compute_distance_gradient",
            &EdgeVertexConstraint::compute_distance_gradient, "", py::arg("V"),
            py::arg("E"), py::arg("F"))
        .def(
            "compute_distance_hessian",
            &EdgeVertexConstraint::compute_distance_hessian, "", py::arg("V"),
            py::arg("E"), py::arg("F"))
        .def(
            "compute_potential", &EdgeVertexConstraint::compute_potential, "",
            py::arg("V"), py::arg("E"), py::arg("F"), py::arg("dhat"))
        .def(
            "compute_potential_gradient",
            &EdgeVertexConstraint::compute_potential_gradient, "", py::arg("V"),
            py::arg("E"), py::arg("F"), py::arg("dhat"))
        .def(
            "compute_potential_hessian",
            &EdgeVertexConstraint::compute_potential_hessian, "", py::arg("V"),
            py::arg("E"), py::arg("F"), py::arg("dhat"),
            py::arg("project_hessian_to_psd"))
        .def_readwrite("multiplicity", &EdgeVertexConstraint::multiplicity, "");

    py::class_<EdgeEdgeConstraint, EdgeEdgeCandidate, CollisionConstraint>(
        m, "EdgeEdgeConstraint")
        .def(
            py::init<long, long, double>(), "", py::arg("edge0_index"),
            py::arg("edge1_index"), py::arg("eps_x"))
        .def(
            py::init<const EdgeEdgeCandidate&, double>(), "",
            py::arg("candidate"), py::arg("eps_x"))
        .def("num_vertices", &EdgeEdgeConstraint::num_vertices, "")
        .def(
            "vertex_indices", &EdgeEdgeConstraint::vertex_indices, "",
            py::arg("E"), py::arg("F"))
        .def(
            "compute_distance", &EdgeEdgeConstraint::compute_distance, "",
            py::arg("V"), py::arg("E"), py::arg("F"))
        .def(
            "compute_distance_gradient",
            &EdgeEdgeConstraint::compute_distance_gradient, "", py::arg("V"),
            py::arg("E"), py::arg("F"))
        .def(
            "compute_distance_hessian",
            &EdgeEdgeConstraint::compute_distance_hessian, "", py::arg("V"),
            py::arg("E"), py::arg("F"))
        .def(
            "compute_potential", &EdgeEdgeConstraint::compute_potential, "",
            py::arg("V"), py::arg("E"), py::arg("F"), py::arg("dhat"))
        .def(
            "compute_potential_gradient",
            &EdgeEdgeConstraint::compute_potential_gradient, "", py::arg("V"),
            py::arg("E"), py::arg("F"), py::arg("dhat"))
        .def(
            "compute_potential_hessian",
            &EdgeEdgeConstraint::compute_potential_hessian, "", py::arg("V"),
            py::arg("E"), py::arg("F"), py::arg("dhat"),
            py::arg("project_hessian_to_psd"))
        .def_readwrite("eps_x", &EdgeEdgeConstraint::eps_x, "");

    py::class_<FaceVertexConstraint, FaceVertexCandidate, CollisionConstraint>(
        m, "FaceVertexConstraint")
        .def(
            py::init<long, long>(), "", py::arg("face_index"),
            py::arg("vertex_index"))
        .def(py::init<const FaceVertexCandidate&>(), "", py::arg("candidate"))
        .def("num_vertices", &FaceVertexConstraint::num_vertices, "")
        .def(
            "vertex_indices", &FaceVertexConstraint::vertex_indices, "",
            py::arg("E"), py::arg("F"))
        .def(
            "compute_distance", &FaceVertexConstraint::compute_distance, "",
            py::arg("V"), py::arg("E"), py::arg("F"))
        .def(
            "compute_distance_gradient",
            &FaceVertexConstraint::compute_distance_gradient, "", py::arg("V"),
            py::arg("E"), py::arg("F"))
        .def(
            "compute_distance_hessian",
            &FaceVertexConstraint::compute_distance_hessian, "", py::arg("V"),
            py::arg("E"), py::arg("F"));

    py::class_<PlaneVertexConstraint, CollisionConstraint>(
        m, "PlaneVertexConstraint")
        .def(
            py::init<const VectorMax3d&, const VectorMax3d&, const long>(), "",
            py::arg("plane_origin"), py::arg("plane_normal"),
            py::arg("vertex_index"))
        .def("num_vertices", &PlaneVertexConstraint::num_vertices, "")
        .def(
            "vertex_indices", &PlaneVertexConstraint::vertex_indices, "",
            py::arg("E"), py::arg("F"))
        .def(
            "compute_distance", &PlaneVertexConstraint::compute_distance, "",
            py::arg("V"), py::arg("E"), py::arg("F"))
        .def(
            "compute_distance_gradient",
            &PlaneVertexConstraint::compute_distance_gradient, "", py::arg("V"),
            py::arg("E"), py::arg("F"))
        .def(
            "compute_distance_hessian",
            &PlaneVertexConstraint::compute_distance_hessian, "", py::arg("V"),
            py::arg("E"), py::arg("F"))
        .def_readwrite("plane_origin", &PlaneVertexConstraint::plane_origin, "")
        .def_readwrite("plane_normal", &PlaneVertexConstraint::plane_normal, "")
        .def_readwrite(
            "vertex_index", &PlaneVertexConstraint::vertex_index, "");

    py::class_<Constraints>(m, "Constraints")
        .def("__len__", &Constraints::size, "")
        .def("num_constraints", &Constraints::num_constraints, "")
        .def("empty", &Constraints::empty, "")
        .def("clear", &Constraints::clear, "")
        .def(
            "__getitem__",
            [](Constraints& self, size_t idx) -> CollisionConstraint* {
                return &self[idx];
            })
        .def_readwrite("vv_constraints", &Constraints::vv_constraints, "")
        .def_readwrite("ev_constraints", &Constraints::ev_constraints, "")
        .def_readwrite("ee_constraints", &Constraints::ee_constraints, "")
        .def_readwrite("fv_constraints", &Constraints::fv_constraints, "")
        .def_readwrite("pv_constraints", &Constraints::pv_constraints, "");
}
