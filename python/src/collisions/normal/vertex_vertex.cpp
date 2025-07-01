#include <common.hpp>

#include <ipc/collisions/normal/vertex_vertex.hpp>

using namespace ipc;

void define_vertex_vertex_normal_collision(py::module_& m)
{
    py::class_<
        VertexVertexNormalCollision, VertexVertexCandidate, NormalCollision>(
        m, "VertexVertexNormalCollision")
        .def(py::init<index_t, index_t>(), "", "vertex0_id"_a, "vertex1_id"_a)
        .def(py::init<const VertexVertexCandidate&>(), "vv_candidate"_a);
    // .def(
    //     py::init<
    //         const index_t, const index_t, const double,
    //         const Eigen::SparseVector<double>&>(),
    //     "vertex0_id"_a, "vertex1_id"_a, "weight"_a,
    //     "weight_gradient"_a);
}
