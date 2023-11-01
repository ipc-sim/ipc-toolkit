#include <common.hpp>

#include <ipc/collisions/edge_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_vertex_constraint(py::module_& m)
{
    py::class_<EdgeVertexConstraint, EdgeVertexCandidate, CollisionConstraint>(
        m, "EdgeVertexConstraint")
        .def(py::init<long, long>(), py::arg("edge_id"), py::arg("vertex_id"))
        .def(py::init<const EdgeVertexCandidate&>(), py::arg("candidate"));
    // .def(
    //     py::init<
    //         const long, const long, const double,
    //         const Eigen::SparseVector<double>&>(),
    //     py::arg("edge_id"), py::arg("vertex_id"), py::arg("weight"),
    //     py::arg("weight_gradient"));
}
