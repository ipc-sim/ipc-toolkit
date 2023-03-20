#include <common.hpp>

#include <ipc/broad_phase/hash_grid.hpp>

namespace py = pybind11;
using namespace ipc;

void define_hash_grid(py::module_& m)
{
    py::class_<HashItem>(m, "HashItem")
        .def(
            py::init<int, int>(),
            "Construct a hash item as a (key, value) pair.", py::arg("key"),
            py::arg("id"))
        .def(
            "__lt__", &HashItem::operator<,
            "Compare HashItems by their keys for sorting.", py::arg("other"))
        .def_readwrite("key", &HashItem::key, "The key of the item.")
        .def_readwrite("id", &HashItem::id, "The value of the item.");

    py::class_<HashGrid, BroadPhase>(m, "HashGrid")
        .def(py::init())
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, double>(&HashGrid::build),
            R"ipc_Qu8mg5v7(
            Build the broad phase for static collision detection.

            Parameters:
                vertices_t0: Vertex positions
                edges: Collision mesh edges
                faces: Collision mesh faces
                inflation_radius: Radius of inflation around all elements.
            )ipc_Qu8mg5v7",
            py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("inflation_radius") = 0)
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, double>(
                &HashGrid::build),
            R"ipc_Qu8mg5v7(
            Build the broad phase for continuous collision detection.

            Parameters:
                vertices_t0: Starting vertices of the vertices.
                vertices_t1: Ending vertices of the vertices.
                edges: Collision mesh edges
                faces: Collision mesh faces
                inflation_radius: Radius of inflation around all elements.
            )ipc_Qu8mg5v7",
            py::arg("vertices_t0"), py::arg("vertices_t1"), py::arg("edges"),
            py::arg("faces"), py::arg("inflation_radius") = 0)
        .def("clear", &HashGrid::clear, "Clear the hash grid.")
        .def(
            "detect_edge_vertex_candidates",
            [](HashGrid& self) {
                std::vector<EdgeVertexCandidate> candidates;
                self.detect_edge_vertex_candidates(candidates);
                return candidates;
            },
            R"ipc_Qu8mg5v7(
            Find the candidate edge-vertex collisisons.

            Returns:
                The candidate edge-vertex collisisons.
            )ipc_Qu8mg5v7")
        .def(
            "detect_edge_edge_candidates",
            [](HashGrid& self) {
                std::vector<EdgeEdgeCandidate> candidates;
                self.detect_edge_edge_candidates(candidates);
                return candidates;
            },
            R"ipc_Qu8mg5v7(
            Find the candidate edge-edge collisions.

            Returns:
                The candidate edge-edge collisisons.
            )ipc_Qu8mg5v7")
        .def(
            "detect_face_vertex_candidates",
            [](HashGrid& self) {
                std::vector<FaceVertexCandidate> candidates;
                self.detect_face_vertex_candidates(candidates);
                return candidates;
            },
            R"ipc_Qu8mg5v7(
            Find the candidate face-vertex collisions.

            Returns:
                The candidate face-vertex collisisons.
            )ipc_Qu8mg5v7")
        .def(
            "detect_edge_face_candidates",
            [](HashGrid& self) {
                std::vector<EdgeFaceCandidate> candidates;
                self.detect_edge_face_candidates(candidates);
                return candidates;
            },
            R"ipc_Qu8mg5v7(
            Find the candidate edge-face intersections.

            Returns:
                The candidate edge-face intersections.
            )ipc_Qu8mg5v7")
        .def("cellSize", &HashGrid::cellSize, "")
        .def("gridSize", &HashGrid::gridSize, "")
        .def("domainMin", &HashGrid::domainMin, "")
        .def("domainMax", &HashGrid::domainMax, "");
}
