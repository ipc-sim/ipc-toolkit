#include <common.hpp>

#include <ipc/collisions/normal/edge_edge.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_edge_normal_collision(py::module_& m)
{
    py::class_<EdgeEdgeNormalCollision, EdgeEdgeCandidate, NormalCollision>(
        m, "EdgeEdgeNormalCollision")
        .def(
            py::init<
                const index_t, const index_t, const double,
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
        //         const index_t, const index_t, const double, const double,
        //         const Eigen::SparseVector<double>&,
        //         const EdgeEdgeDistanceType>(),
        //     py::arg("edge0_id"), py::arg("edge1_id"), py::arg("eps_x"),
        //     py::arg("weight"), py::arg("weight_gradient"),
        //     py::arg("dtype") = EdgeEdgeDistanceType::AUTO)
        .def("__eq__", &EdgeEdgeNormalCollision::operator==, py::arg("other"))
        .def("__ne__", &EdgeEdgeNormalCollision::operator!=, py::arg("other"))
        .def("__lt__", &EdgeEdgeNormalCollision::operator<, py::arg("other"))
        .def_readwrite(
            "eps_x", &EdgeEdgeNormalCollision::eps_x,
            "Mollifier activation threshold.")
        .def_readwrite(
            "dtype", &EdgeEdgeNormalCollision::dtype,
            R"ipc_Qu8mg5v7(
            Cached distance type.

            Some EE collisions are mollified EV or VV collisions.
            )ipc_Qu8mg5v7");
}
