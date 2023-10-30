#include <common.hpp>

#include <ipc/friction/constraints/edge_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_vertex_friction_constraint(py::module_& m)
{
    py::class_<
        EdgeVertexFrictionConstraint, EdgeVertexCandidate, FrictionConstraint>(
        m, "EdgeVertexFrictionConstraint")
        .def(py::init<const EdgeVertexConstraint&>(), py::arg("constraint"))
        .def(
            py::init<
                const EdgeVertexConstraint&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double>(),
            py::arg("constraint"), py::arg("vertices"), py::arg("edges"),
            py::arg("faces"), py::arg("dhat"), py::arg("barrier_stiffness"));
}
