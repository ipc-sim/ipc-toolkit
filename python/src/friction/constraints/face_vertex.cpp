#include <common.hpp>

#include <ipc/friction/constraints/face_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_face_vertex_friction_constraint(py::module_& m)
{
    py::class_<
        FaceVertexFrictionConstraint, FaceVertexCandidate, FrictionConstraint>(
        m, "FaceVertexFrictionConstraint")
        .def(py::init<const FaceVertexConstraint&>(), py::arg("constraint"))
        .def(
            py::init<
                const FaceVertexConstraint&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double>(),
            py::arg("constraint"), py::arg("vertices"), py::arg("edges"),
            py::arg("faces"), py::arg("dhat"), py::arg("barrier_stiffness"));
}
