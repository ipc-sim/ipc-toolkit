#include <common.hpp>

#include <ipc/collisions/normal/face_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_face_vertex_normal_collision(py::module_& m)
{
    py::class_<FaceVertexNormalCollision, FaceVertexCandidate, NormalCollision>(
        m, "FaceVertexNormalCollision")
        .def(
            py::init<long, long>(), "", py::arg("face_id"),
            py::arg("vertex_id"))
        .def(py::init<const FaceVertexCandidate&>(), py::arg("candidate"));
    // .def(
    //     py::init<
    //         const long, const long, const double,
    //         const Eigen::SparseVector<double>&>(),
    //     py::arg("face_id"), py::arg("vertex_id"), py::arg("weight"),
    //     py::arg("weight_gradient"));
}
