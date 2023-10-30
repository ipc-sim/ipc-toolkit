#include <common.hpp>

#include <ipc/broad_phase/bvh.hpp>

namespace py = pybind11;
using namespace ipc;

void define_bvh(py::module_& m)
{
    py::class_<BVH, BroadPhase>(m, "BVH")
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, const double>(&BVH::build),
            R"ipc_Qu8mg5v7(
            Build the broad phase for static collision detection.

            Parameters:
                V0: Positions of the vertices.
                E: Edges of the mesh.
                F: Faces of the mesh.
                inflation_radius: Radius of inflation around all elements.
            )ipc_Qu8mg5v7",
            py::arg("V"), py::arg("E"), py::arg("F"),
            py::arg("inflation_radius") = 0)
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double>(
                &BVH::build),
            R"ipc_Qu8mg5v7(
            Build the broad phase for continuous collision detection.

            Parameters:
                V0: Starting positions of the vertices.
                V1: Ending positions of the vertices.
                E: Edges of the mesh.
                F: Faces of the mesh.
                inflation_radius: Radius of inflation around all elements.
            )ipc_Qu8mg5v7",
            py::arg("V0"), py::arg("V1"), py::arg("E"), py::arg("F"),
            py::arg("inflation_radius") = 0)
        .def("clear", &BVH::clear, "Clear any built data.")
        .def(
            "detect_vertex_vertex_candidates",
            [](BVH& self) {
                std::vector<VertexVertexCandidate> candidates;
                self.detect_vertex_vertex_candidates(candidates);
                return candidates;
            },
            R"ipc_Qu8mg5v7(
            Find the candidate vertex-vertex collisions.

            Returns:
                The candidate vertex-vertex collisions.
            )ipc_Qu8mg5v7")
        .def(
            "detect_edge_vertex_candidates",
            [](BVH& self) {
                std::vector<EdgeVertexCandidate> candidates;
                self.detect_edge_vertex_candidates(candidates);
                return candidates;
            },
            R"ipc_Qu8mg5v7(
            Find the candidate edge-vertex collisions.

            Returns:
                The candidate edge-vertex collisions.
            )ipc_Qu8mg5v7")
        .def(
            "detect_edge_edge_candidates",
            [](BVH& self) {
                std::vector<EdgeEdgeCandidate> candidates;
                self.detect_edge_edge_candidates(candidates);
                return candidates;
            },
            R"ipc_Qu8mg5v7(
            Find the candidate edge-edge collisions.

            Returns:
                The candidate edge-edge collisions.
            )ipc_Qu8mg5v7")
        .def(
            "detect_face_vertex_candidates",
            [](BVH& self) {
                std::vector<FaceVertexCandidate> candidates;
                self.detect_face_vertex_candidates(candidates);
                return candidates;
            },
            R"ipc_Qu8mg5v7(
            Find the candidate face-vertex collisions.

            Returns:
                The candidate face-vertex collisions.
            )ipc_Qu8mg5v7")
        .def(
            "detect_edge_face_candidates",
            [](BVH& self) {
                std::vector<EdgeFaceCandidate> candidates;
                self.detect_edge_face_candidates(candidates);
                return candidates;
            },
            R"ipc_Qu8mg5v7(
            Find the candidate edge-face intersections.

            Returns:
                The candidate edge-face intersections.
            )ipc_Qu8mg5v7");
}
