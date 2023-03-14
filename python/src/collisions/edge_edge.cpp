#include <common.hpp>

#include <ipc/collisions/edge_edge.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_edge_constraint(py::module_& m)
{
    py::class_<EdgeEdgeConstraint, EdgeEdgeCandidate, CollisionConstraint>(
        m, "EdgeEdgeConstraint")
        .def(
            py::init<long, long, double>(), "", py::arg("edge0_id"),
            py::arg("edge1_id"), py::arg("eps_x"))
        .def(
            py::init<const EdgeEdgeCandidate&, double>(), "",
            py::arg("candidate"), py::arg("eps_x"))
        .def(
            "compute_potential", &EdgeEdgeConstraint::compute_potential, "",
            py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("dhat"))
        .def(
            "compute_potential_gradient",
            &EdgeEdgeConstraint::compute_potential_gradient, "",
            py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("dhat"))
        .def(
            "compute_potential_hessian",
            &EdgeEdgeConstraint::compute_potential_hessian, "",
            py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("dhat"), py::arg("project_hessian_to_psd"))
        .def_readwrite("eps_x", &EdgeEdgeConstraint::eps_x, "");
}
