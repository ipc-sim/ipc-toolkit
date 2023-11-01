#include <common.hpp>

#include <ipc/friction/constraints/vertex_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_vertex_vertex_friction_constraint(py::module_& m)
{
    py::class_<
        VertexVertexFrictionConstraint, VertexVertexCandidate,
        FrictionConstraint>(m, "VertexVertexFrictionConstraint")
        .def(py::init<const VertexVertexConstraint&>(), py::arg("constraint"))
        .def(
            py::init<
                const VertexVertexConstraint&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double>(),
            py::arg("constraint"), py::arg("vertices"), py::arg("edges"),
            py::arg("faces"), py::arg("dhat"), py::arg("barrier_stiffness"));
}
