#include <common.hpp>

#include <ipc/broad_phase/spatial_hash.hpp>

using namespace ipc;

void define_spatial_hash(py::module_& m)
{
    py::class_<SpatialHash, BroadPhase, std::shared_ptr<SpatialHash>>(
        m, "SpatialHash")
        .def(py::init())
        .def(
            py::init<
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>, double, double>(),
            "vertices"_a, "edges"_a, "faces"_a, "inflation_radius"_a = 0,
            "voxel_size"_a = -1)
        .def(
            py::init<
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>, double, double>(),
            "vertices_t0"_a, "vertices_t1"_a, "edges"_a, "faces"_a,
            "inflation_radius"_a = 0, "voxel_size"_a = -1)
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
            "inflation_radius"_a = 0, "voxel_size"_a = -1)
        .def("clear", &SpatialHash::clear)
        .def(
            "is_vertex_index", &SpatialHash::is_vertex_index,
            "Check if primitive index refers to a vertex.", "idx"_a)
        .def(
            "is_edge_index", &SpatialHash::is_edge_index,
            "Check if primitive index refers to an edge.", "idx"_a)
        .def(
            "is_triangle_index", &SpatialHash::is_triangle_index,
            "Check if primitive index refers to a triangle.", "idx"_a)
        .def(
            "to_edge_index", &SpatialHash::to_edge_index,
            "Convert a primitive index to an edge index.", "idx"_a)
        .def(
            "to_triangle_index", &SpatialHash::to_triangle_index,
            "Convert a primitive index to a triangle index.", "idx"_a)
        .def_readwrite(
            "left_bottom_corner", &SpatialHash::left_bottom_corner,
            "The left bottom corner of the world bounding box.")
        .def_readwrite(
            "right_top_corner", &SpatialHash::right_top_corner,
            "The right top corner of the world bounding box.")
        .def_readwrite(
            "voxel_count", &SpatialHash::voxel_count,
            "The number of voxels in each dimension.")
        .def_readwrite(
            "one_div_voxelSize", &SpatialHash::one_div_voxelSize,
            "1.0 / voxel_size")
        .def_readwrite(
            "voxel_count_0x1", &SpatialHash::voxel_count_0x1,
            "The number of voxels in the first two dimensions.")
        .def_readwrite("edge_start_ind", &SpatialHash::edge_start_ind)
        .def_readwrite("tri_start_ind", &SpatialHash::tri_start_ind)
        .def_readwrite(
            "voxel_to_primitives", &SpatialHash::voxel_to_primitives,
            "Map from voxel index to the primitive indices it contains.")
        .def_readwrite(
            "point_to_voxels", &SpatialHash::point_to_voxels,
            "Map from point index to the voxel indices it occupies.")
        .def_readwrite(
            "edge_to_voxels", &SpatialHash::edge_to_voxels,
            "Map from edge index to the voxel indices it occupies.")
        .def_readwrite(
            "face_to_voxels", &SpatialHash::face_to_voxels,
            "Map from face index to the voxel indices it occupies.");
}
