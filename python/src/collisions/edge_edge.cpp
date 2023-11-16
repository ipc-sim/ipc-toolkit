#include <common.hpp>

#include <ipc/collisions/edge_edge.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_edge_constraint(py::module_& m)
{
    py::class_<EdgeEdgeConstraint, EdgeEdgeCandidate, CollisionConstraint>(
        m, "EdgeEdgeConstraint")
        .def(
            py::init<
                const long, const long, const double,
                const EdgeEdgeDistanceType>(),
            py::arg("edge0_id"), py::arg("edge1_id"), py::arg("eps_x"),
            py::arg("dtype") = EdgeEdgeDistanceType::AUTO)
        .def(
            py::init<
                const EdgeEdgeCandidate&, const double,
                const EdgeEdgeDistanceType>(),
            py::arg("candidate"), py::arg("eps_x"),
            py::arg("dtype") = EdgeEdgeDistanceType::AUTO)
        // .def(
        //     py::init<
        //         const long, const long, const double, const double,
        //         const Eigen::SparseVector<double>&,
        //         const EdgeEdgeDistanceType>(),
        //     py::arg("edge0_id"), py::arg("edge1_id"), py::arg("eps_x"),
        //     py::arg("weight"), py::arg("weight_gradient"),
        //     py::arg("dtype") = EdgeEdgeDistanceType::AUTO)
        .def("__eq__", &EdgeEdgeConstraint::operator==, py::arg("other"))
        .def("__ne__", &EdgeEdgeConstraint::operator!=, py::arg("other"))
        .def("__lt__", &EdgeEdgeConstraint::operator<, py::arg("other"))
        .def_readwrite(
            "eps_x", &EdgeEdgeConstraint::eps_x,
            "Mollifier activation threshold.")
        .def_readwrite(
            "dtype", &EdgeEdgeConstraint::dtype,
            R"ipc_Qu8mg5v7(
            Cached distance type.

            Some EE constraints are mollified EV or VV constraints.
            )ipc_Qu8mg5v7");
}
