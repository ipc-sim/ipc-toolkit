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
            "", py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("inflation_radius") = 0, py::arg("voxel_size") = -1)
        .def(
            py::init<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, double,
                double>(),
            "", py::arg("vertices_t0"), py::arg("vertices_t1"),
            py::arg("edges"), py::arg("faces"), py::arg("inflation_radius") = 0,
            py::arg("voxel_size") = -1)
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, double, double>(&SpatialHash::build),
            "", py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("inflation_radius") = 0, py::arg("voxel_size") = -1)
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, double, double>(
                &SpatialHash::build),
            "", py::arg("vertices_t0"), py::arg("vertices_t1"),
            py::arg("edges"), py::arg("faces"), py::arg("inflation_radius") = 0,
            py::arg("voxel_size") = -1)
        .def("clear", &SpatialHash::clear)
        .def(
            "detect_edge_vertex_candidates",
            [](SpatialHash& self) {
                std::vector<EdgeVertexCandidate> candidates;
                self.detect_edge_vertex_candidates(candidates);
                return candidates;
            },
            "Find the candidate edge-vertex collisisons.")
        .def(
            "detect_edge_edge_candidates",
            [](SpatialHash& self) {
                std::vector<EdgeEdgeCandidate> candidates;
                self.detect_edge_edge_candidates(candidates);
                return candidates;
            },
            "Find the candidate edge-edge collisions.")
        .def(
            "detect_face_vertex_candidates",
            [](SpatialHash& self) {
                std::vector<FaceVertexCandidate> candidates;
                self.detect_face_vertex_candidates(candidates);
                return candidates;
            },
            "Find the candidate face-vertex collisions.")
        .def(
            "detect_edge_face_candidates",
            [](SpatialHash& self) {
                std::vector<EdgeFaceCandidate> candidates;
                self.detect_edge_face_candidates(candidates);
                return candidates;
            },
            "Find the candidate edge-face intersections.")
        .def_readwrite("left_bottom_corner", &SpatialHash::left_bottom_corner)
        .def_readwrite("right_top_corner", &SpatialHash::right_top_corner)
        .def_readwrite("voxel_count", &SpatialHash::voxel_count)
        .def_readwrite("one_div_voxelSize", &SpatialHash::one_div_voxelSize)
        .def_readwrite("voxel_count_0x1", &SpatialHash::voxel_count_0x1)
        .def_readwrite("edge_start_ind", &SpatialHash::edge_start_ind)
        .def_readwrite("tri_start_ind", &SpatialHash::tri_start_ind)
        .def_readwrite("voxel", &SpatialHash::voxel)
        .def_readwrite(
            "point_and_edge_occupancy", &SpatialHash::point_and_edge_occupancy);
}
