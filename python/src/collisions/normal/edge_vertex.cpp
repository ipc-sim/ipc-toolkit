#include <common.hpp>

#include <ipc/collisions/normal/edge_vertex.hpp>

using namespace ipc;

void define_edge_vertex_normal_collision(py::module_& m)
{
    py::class_<EdgeVertexNormalCollision, EdgeVertexCandidate, NormalCollision>(
        m, "EdgeVertexNormalCollision")
        .def(py::init<index_t, index_t>(), "edge_id"_a, "vertex_id"_a)
        .def(py::init<const EdgeVertexCandidate&>(), "candidate"_a);
    // .def(
    //     py::init<
    //         const index_t, const index_t, const double,
    //         const Eigen::SparseVector<double>&>(),
    //     "edge_id"_a, "vertex_id"_a, "weight"_a,
    //     "weight_gradient"_a);
}
