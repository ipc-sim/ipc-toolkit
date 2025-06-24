#include <common.hpp>

#include <ipc/collisions/normal/face_vertex.hpp>

using namespace ipc;

void define_face_vertex_normal_collision(py::module_& m)
{
    py::class_<FaceVertexNormalCollision, FaceVertexCandidate, NormalCollision>(
        m, "FaceVertexNormalCollision")
        .def(py::init<index_t, index_t>(), "", "face_id"_a, "vertex_id"_a)
        .def(py::init<const FaceVertexCandidate&>(), "candidate"_a);
    // .def(
    //     py::init<
    //         const index_t, const index_t, const double,
    //         const Eigen::SparseVector<double>&>(),
    //     "face_id"_a, "vertex_id"_a, "weight"_a,
    //     "weight_gradient"_a);
}
