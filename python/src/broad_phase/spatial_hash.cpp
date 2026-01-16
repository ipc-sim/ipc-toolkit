#include <common.hpp>

#include <ipc/broad_phase/spatial_hash.hpp>

using namespace ipc;

void define_spatial_hash(py::module_& m)
{
    py::class_<SpatialHash, BroadPhase, std::shared_ptr<SpatialHash>>(
        m, "SpatialHash")
        .def(py::init())
        .def(
            "build",
            py::overload_cast<
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>, double, double>(
                &SpatialHash::build),
            "vertices"_a, "edges"_a, "faces"_a, "inflation_radius"_a = 0,
            "voxel_size"_a = -1)
        .def(
            "build",
            py::overload_cast<
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>, double, double>(
                &SpatialHash::build),
            "vertices_t0"_a, "vertices_t1"_a, "edges"_a, "faces"_a,
            "inflation_radius"_a = 0, "voxel_size"_a = -1);
}
