#include <common.hpp>

#include <ipc/friction/collisions/edge_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_vertex_friction_collision(py::module_& m)
{
    py::class_<
        EdgeVertexFrictionCollision, EdgeVertexCandidate, FrictionCollision>(
        m, "EdgeVertexFrictionCollision")
        .def(py::init<const EdgeVertexCollision&>(), py::arg("collision"))
        .def(
            py::init<
                const EdgeVertexCollision&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double>(),
            py::arg("collision"), py::arg("vertices"), py::arg("edges"),
            py::arg("faces"), py::arg("dhat"), py::arg("barrier_stiffness"));
}
