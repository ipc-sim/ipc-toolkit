#include <common.hpp>

#include <ipc/friction/collisions/edge_edge.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_edge_friction_collision(py::module_& m)
{
    py::class_<EdgeEdgeFrictionCollision, EdgeEdgeCandidate, FrictionCollision>(
        m, "EdgeEdgeFrictionCollision")
        .def(py::init<const EdgeEdgeCollision&>(), py::arg("collision"))
        .def(
            py::init<
                const EdgeEdgeCollision&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double>(),
            py::arg("collision"), py::arg("vertices"), py::arg("edges"),
            py::arg("faces"), py::arg("dhat"), py::arg("barrier_stiffness"));
}
