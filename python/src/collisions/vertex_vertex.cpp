#include <common.hpp>

#include <ipc/collisions/vertex_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_vertex_vertex_collision(py::module_& m)
{
    py::class_<VertexVertexCollision, VertexVertexCandidate, Collision>(
        m, "VertexVertexCollision")
        .def(
            py::init<long, long>(), "", py::arg("vertex0_id"),
            py::arg("vertex1_id"))
        .def(py::init<const VertexVertexCandidate&>(), py::arg("vv_candidate"));
    // .def(
    //     py::init<
    //         const long, const long, const double,
    //         const Eigen::SparseVector<double>&>(),
    //     py::arg("vertex0_id"), py::arg("vertex1_id"), py::arg("weight"),
    //     py::arg("weight_gradient"));
}
