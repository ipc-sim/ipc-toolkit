#include <common.hpp>

#include <ipc/collisions/face_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_face_vertex_collision(py::module_& m)
{
    py::class_<FaceVertexCollision, FaceVertexCandidate, Collision>(
        m, "FaceVertexCollision")
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
