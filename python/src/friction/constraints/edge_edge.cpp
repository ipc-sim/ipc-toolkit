#include <common.hpp>

#include <ipc/friction/constraints/edge_edge.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_edge_friction_constraint(py::module_& m)
{
    py::class_<
        EdgeEdgeFrictionConstraint, EdgeEdgeCandidate, FrictionConstraint>(
        m, "EdgeEdgeFrictionConstraint")
        .def(py::init<const EdgeEdgeConstraint&>(), py::arg("constraint"))
        .def(
            py::init<
                const EdgeEdgeConstraint&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double>(),
            py::arg("constraint"), py::arg("vertices"), py::arg("edges"),
            py::arg("faces"), py::arg("dhat"), py::arg("barrier_stiffness"));
}
