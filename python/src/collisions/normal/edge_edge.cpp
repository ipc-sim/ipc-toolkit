#include <common.hpp>

#include <ipc/collisions/normal/edge_edge.hpp>

using namespace ipc;

void define_edge_edge_normal_collision(py::module_& m)
{
    py::class_<EdgeEdgeNormalCollision, EdgeEdgeCandidate, NormalCollision>(
        m, "EdgeEdgeNormalCollision")
        .def(
            py::init<
                const index_t, const index_t, const double,
                const EdgeEdgeDistanceType>(),
            "edge0_id"_a, "edge1_id"_a, "eps_x"_a,
            "dtype"_a = EdgeEdgeDistanceType::AUTO)
        .def(
            py::init<
                const EdgeEdgeCandidate&, const double,
                const EdgeEdgeDistanceType>(),
            "candidate"_a, "eps_x"_a, "dtype"_a = EdgeEdgeDistanceType::AUTO)
        // .def(
        //     py::init<
        //         const index_t, const index_t, const double, const double,
        //         const Eigen::SparseVector<double>&,
        //         const EdgeEdgeDistanceType>(),
        //     "edge0_id"_a, "edge1_id"_a, "eps_x"_a,
        //     "weight"_a, "weight_gradient"_a,
        //     "dtype"_a = EdgeEdgeDistanceType::AUTO)
        .def("__eq__", &EdgeEdgeNormalCollision::operator==, "other"_a)
        .def("__ne__", &EdgeEdgeNormalCollision::operator!=, "other"_a)
        .def("__lt__", &EdgeEdgeNormalCollision::operator<, "other"_a)
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
