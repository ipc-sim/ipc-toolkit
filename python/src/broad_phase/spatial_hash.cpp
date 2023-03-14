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
            py::arg("inflation_radius") = 0, py::arg("voxelSize") = -1)
        .def(
            py::init<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, double,
                double>(),
            "", py::arg("vertices_t0"), py::arg("vertices_t1"),
            py::arg("edges"), py::arg("faces"), py::arg("inflation_radius") = 0,
            py::arg("voxelSize") = -1)
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, double>(&SpatialHash::build),
            "", py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("inflation_radius") = 0)
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, double>(
                &SpatialHash::build),
            "", py::arg("vertices_t0"), py::arg("vertices_t1"),
            py::arg("edges"), py::arg("faces"), py::arg("inflation_radius") = 0)
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, double, double>(&SpatialHash::build),
            "", py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("inflation_radius"), py::arg("voxelSize"))
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, double, double>(
                &SpatialHash::build),
            "", py::arg("vertices_t0"), py::arg("vertices_t1"),
            py::arg("edges"), py::arg("faces"), py::arg("inflation_radius"),
            py::arg("voxelSize"))
        .def("clear", &SpatialHash::clear, "")
        .def(
            "queryPointForTriangles",
            [](SpatialHash& self, const VectorMax3d& p, double radius = 0) {
                unordered_set<int> triInds;
                self.queryPointForTriangles(p, triInds, radius);
                return triInds;
            },
            "", py::arg("p"), py::arg("radius") = 0)
        .def(
            "queryPointForTriangles",
            [](SpatialHash& self, const VectorMax3d& p_t0,
               const VectorMax3d& p_t1, double radius = 0) {
                unordered_set<int> triInds;
                self.queryPointForTriangles(p_t0, p_t1, triInds, radius);
                return triInds;
            },
            "", py::arg("p_t0"), py::arg("p_t1"), py::arg("radius") = 0)
        .def(
            "queryPointForPrimitives",
            [](SpatialHash& self, const VectorMax3d& p_t0,
               const VectorMax3d& p_t1, double radius = 0) {
                unordered_set<int> vertInds;
                unordered_set<int> edgeInds;
                unordered_set<int> triInds;
                self.queryPointForPrimitives(
                    p_t0, p_t1, vertInds, edgeInds, triInds, radius);
                return std::make_tuple(vertInds, edgeInds, triInds);
            },
            "", py::arg("p_t0"), py::arg("p_t1"), py::arg("radius") = 0)
        .def(
            "queryEdgeForPE",
            [](SpatialHash& self, const VectorMax3d& e0, const VectorMax3d& e1,
               double radius = 0) {
                std::vector<int> vertInds;
                std::vector<int> edgeInds;
                self.queryEdgeForPE(e0, e1, vertInds, edgeInds, radius);
                return std::make_tuple(vertInds, edgeInds);
            },
            "", py::arg("e0"), py::arg("e1"), py::arg("radius") = 0)
        .def(
            "queryEdgeForEdges",
            [](SpatialHash& self, const VectorMax3d& ea0,
               const VectorMax3d& ea1, double radius = 0, int eai = -1) {
                std::vector<int> edgeInds;
                self.queryEdgeForEdges(ea0, ea1, edgeInds, radius, eai);
                return edgeInds;
            },
            "", py::arg("ea0"), py::arg("ea1"), py::arg("radius") = 0,
            py::arg("eai") = -1)
        .def(
            "queryEdgeForEdgesWithBBoxCheck",
            [](SpatialHash& self, const Eigen::MatrixXd& vertices,
               const Eigen::MatrixXi& edges, const VectorMax3d& ea0,
               const VectorMax3d& ea1, double radius = 0, int eai = -1) {
                std::vector<int> edgeInds;
                self.queryEdgeForEdgesWithBBoxCheck(
                    vertices, edges, ea0, ea1, edgeInds, radius, eai);
                return edgeInds;
            },
            "", py::arg("vertices"), py::arg("edges"), py::arg("ea0"),
            py::arg("ea1"), py::arg("radius") = 0, py::arg("eai") = -1)
        .def(
            "queryEdgeForEdges",
            [](SpatialHash& self, const VectorMax3d& ea0_t0,
               const VectorMax3d& ea1_t0, const VectorMax3d& ea0_t1,
               const VectorMax3d& ea1_t1, double radius = 0, int eai = -1) {
                std::vector<int> edgeInds;
                self.queryEdgeForEdges(
                    ea0_t0, ea1_t0, ea0_t1, ea1_t1, edgeInds, radius, eai);
                return edgeInds;
            },
            "", py::arg("ea0_t0"), py::arg("ea1_t0"), py::arg("ea0_t1"),
            py::arg("ea1_t1"), py::arg("radius") = 0, py::arg("eai") = -1)
        .def(
            "queryTriangleForPoints",
            [](SpatialHash& self, const VectorMax3d& t0, const VectorMax3d& t1,
               const VectorMax3d& t2, double radius = 0) {
                unordered_set<int> pointInds;
                self.queryTriangleForPoints(t0, t1, t2, pointInds, radius);
                return pointInds;
            },
            "", py::arg("t0"), py::arg("t1"), py::arg("t2"),
            py::arg("radius") = 0)
        .def(
            "queryTriangleForPoints",
            [](SpatialHash& self, const VectorMax3d& t0_t0,
               const VectorMax3d& t1_t0, const VectorMax3d& t2_t0,
               const VectorMax3d& t0_t1, const VectorMax3d& t1_t1,
               const VectorMax3d& t2_t1, double radius = 0) {
                unordered_set<int> pointInds;
                self.queryTriangleForPoints(
                    t0_t0, t1_t0, t2_t0, t0_t1, t1_t1, t2_t1, pointInds,
                    radius);
                return pointInds;
            },
            "", py::arg("t0_t0"), py::arg("t1_t0"), py::arg("t2_t0"),
            py::arg("t0_t1"), py::arg("t1_t1"), py::arg("t2_t1"),
            py::arg("radius") = 0)
        .def(
            "queryTriangleForEdges",
            [](SpatialHash& self, const VectorMax3d& t0, const VectorMax3d& t1,
               const VectorMax3d& t2, double radius = 0) {
                unordered_set<int> edgeInds;
                self.queryTriangleForEdges(t0, t1, t2, edgeInds, radius);
                return edgeInds;
            },
            "", py::arg("t0"), py::arg("t1"), py::arg("t2"),
            py::arg("radius") = 0)
        .def(
            "queryEdgeForTriangles",
            [](SpatialHash& self, const VectorMax3d& e0, const VectorMax3d& e1,
               double radius = 0) {
                unordered_set<int> triInds;
                self.queryEdgeForTriangles(e0, e1, triInds, radius);
                return triInds;
            },
            "", py::arg("e0"), py::arg("e1"), py::arg("radius") = 0)
        .def(
            "queryPointForPrimitives",
            [](SpatialHash& self, int vi) {
                unordered_set<int> vertInds;
                unordered_set<int> edgeInds;
                unordered_set<int> triInds;
                self.queryPointForPrimitives(vi, vertInds, edgeInds, triInds);
                return std::make_tuple(vertInds, edgeInds, triInds);
            },
            "", py::arg("vi"))
        .def(
            "queryPointForEdges",
            [](SpatialHash& self, int vi) {
                unordered_set<int> edgeInds;
                self.queryPointForEdges(vi, edgeInds);
                return edgeInds;
            },
            "", py::arg("vi"))
        .def(
            "queryPointForTriangles",
            [](SpatialHash& self, int vi) {
                unordered_set<int> triInds;
                self.queryPointForTriangles(vi, triInds);
                return triInds;
            },
            "", py::arg("vi"))
        .def(
            "queryEdgeForEdges",
            [](SpatialHash& self, int eai) {
                unordered_set<int> edgeInds;
                self.queryEdgeForEdges(eai, edgeInds);
                return edgeInds;
            },
            "", py::arg("eai"))
        .def(
            "queryEdgeForEdgesWithBBoxCheck",
            [](SpatialHash& self, const Eigen::MatrixXd& vertices_t0,
               const Eigen::MatrixXd& vertices_t1, const Eigen::MatrixXi& edges,
               int eai) {
                unordered_set<int> edgeInds;
                self.queryEdgeForEdgesWithBBoxCheck(
                    vertices_t0, vertices_t1, edges, eai, edgeInds);
                return edgeInds;
            },
            "", py::arg("vertices_t0"), py::arg("vertices_t1"),
            py::arg("edges"), py::arg("eai"))
        .def(
            "queryEdgeForTriangles",
            [](SpatialHash& self, int ei) {
                unordered_set<int> triInds;
                self.queryEdgeForTriangles(ei, triInds);
                return triInds;
            },
            "", py::arg("ei"))
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
        .def_readwrite("leftBottomCorner", &SpatialHash::leftBottomCorner, "")
        .def_readwrite("rightTopCorner", &SpatialHash::rightTopCorner, "")
        .def_readwrite("voxelCount", &SpatialHash::voxelCount, "")
        .def_readwrite("one_div_voxelSize", &SpatialHash::one_div_voxelSize, "")
        .def_readwrite("voxelCount0x1", &SpatialHash::voxelCount0x1, "")
        .def_readwrite("edgeStartInd", &SpatialHash::edgeStartInd, "")
        .def_readwrite("triStartInd", &SpatialHash::triStartInd, "")
        .def_readwrite("voxel", &SpatialHash::voxel, "")
        .def_readwrite(
            "pointAndEdgeOccupancy", &SpatialHash::pointAndEdgeOccupancy, "");
}
