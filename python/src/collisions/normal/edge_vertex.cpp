#include <common.hpp>

#include <ipc/collisions/normal/edge_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_vertex_normal_collision(py::module_& m)
{
    py::class_<EdgeVertexNormalCollision, EdgeVertexCandidate, NormalCollision>(
        m, "EdgeVertexNormalCollision")
        .def(
            py::init<index_t, index_t>(), py::arg("edge_id"),
            py::arg("vertex_id"))
        .def(py::init<const EdgeVertexCandidate&>(), py::arg("candidate"));
    // .def(
    //     py::init<
    //         const index_t, const index_t, const double,
    //         const Eigen::SparseVector<double>&>(),
    //     py::arg("edge_id"), py::arg("vertex_id"), py::arg("weight"),
    //     py::arg("weight_gradient"));
}
