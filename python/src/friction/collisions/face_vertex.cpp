#include <common.hpp>

#include <ipc/friction/collisions/face_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_face_vertex_friction_collision(py::module_& m)
{
    py::class_<
        FaceVertexFrictionCollision, FaceVertexCandidate, FrictionCollision>(
        m, "FaceVertexFrictionCollision")
        .def(py::init<const FaceVertexCollision&>(), py::arg("collision"))
        .def(
            py::init<
                const FaceVertexCollision&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double>(),
            py::arg("collision"), py::arg("vertices"), py::arg("edges"),
            py::arg("faces"), py::arg("dhat"), py::arg("barrier_stiffness"));
}
