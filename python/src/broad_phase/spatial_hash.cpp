#include <common.hpp>

#include <ipc/broad_phase/spatial_hash.hpp>

namespace py = pybind11;
using namespace ipc;

void define_spatial_hash(py::module_& m)
{
    py::class_<SpatialHash, BroadPhase>(m, "SpatialHash")
        .def(py::init())
        .def(
            py::init<
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, double, double>(),
            py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("inflation_radius") = 0, py::arg("voxel_size") = -1)
        .def(
            py::init<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, double,
                double>(),
            py::arg("vertices_t0"), py::arg("vertices_t1"), py::arg("edges"),
            py::arg("faces"), py::arg("inflation_radius") = 0,
            py::arg("voxel_size") = -1)
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, double, double>(&SpatialHash::build),
            py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("inflation_radius") = 0, py::arg("voxel_size") = -1)
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, double, double>(
                &SpatialHash::build),
            py::arg("vertices_t0"), py::arg("vertices_t1"), py::arg("edges"),
            py::arg("faces"), py::arg("inflation_radius") = 0,
            py::arg("voxel_size") = -1)
        .def("clear", &SpatialHash::clear)
        .def(
            "is_vertex_index", &SpatialHash::is_vertex_index,
            "Check if primitive index refers to a vertex.", py::arg("idx"))
        .def(
            "is_edge_index", &SpatialHash::is_edge_index,
            "Check if primitive index refers to an edge.", py::arg("idx"))
        .def(
            "is_triangle_index", &SpatialHash::is_triangle_index,
            "Check if primitive index refers to a triangle.", py::arg("idx"))
        .def(
            "to_edge_index", &SpatialHash::to_edge_index,
            "Convert a primitive index to an edge index.", py::arg("idx"))
        .def(
            "to_triangle_index", &SpatialHash::to_triangle_index,
            "Convert a primitive index to a triangle index.", py::arg("idx"))
        .def_readwrite("left_bottom_corner", &SpatialHash::left_bottom_corner)
        .def_readwrite("right_top_corner", &SpatialHash::right_top_corner)
        .def_readwrite("voxel_count", &SpatialHash::voxel_count)
        .def_readwrite("one_div_voxelSize", &SpatialHash::one_div_voxelSize)
        .def_readwrite("voxel_count_0x1", &SpatialHash::voxel_count_0x1)
        .def_readwrite("edge_start_ind", &SpatialHash::edge_start_ind)
        .def_readwrite("tri_start_ind", &SpatialHash::tri_start_ind)
        .def_readwrite("voxel", &SpatialHash::voxel_to_primitives)
        .def_readwrite("point_to_voxels", &SpatialHash::point_to_voxels)
        .def_readwrite("edge_to_voxels", &SpatialHash::edge_to_voxels)
        .def_readwrite("face_to_voxels", &SpatialHash::face_to_voxels);
}
