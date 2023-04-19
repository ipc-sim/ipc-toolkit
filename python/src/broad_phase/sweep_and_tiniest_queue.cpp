#include <common.hpp>

#include <ipc/broad_phase/sweep_and_tiniest_queue.hpp>

namespace py = pybind11;
using namespace ipc;

void define_sweep(py::module_& m)
{
    py::class_<CopyMeshBroadPhase, BroadPhase>(m, "CopyMeshBroadPhase");

    py::class_<SweepAndTiniestQueue, CopyMeshBroadPhase>(
        m, "SweepAndTiniestQueue")
        .def(py::init())
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, double>(&SweepAndTiniestQueue::build),
            R"ipc_Qu8mg5v7(
            Build the broad phase for static collision detection.

            Parameters:
                vertices: Vertex positions
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
                &SweepAndTiniestQueue::build),
            R"ipc_Qu8mg5v7(
            Build the broad phase for continuous collision detection.

            Parameters:
                vertices_t0: Starting vertex positions
                vertices_t1: Ending vertex positions
                edges: Collision mesh edges
                faces: Collision mesh faces
                inflation_radius: Radius of inflation around all elements.
            )ipc_Qu8mg5v7",
            py::arg("vertices_t0"), py::arg("vertices_t1"), py::arg("edges"),
            py::arg("faces"), py::arg("inflation_radius") = 0)
        .def("clear", &SweepAndTiniestQueue::clear, "Clear any built data.")
        .def(
            "detect_edge_vertex_candidates",
            [](SweepAndTiniestQueue& self) {
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
            [](SweepAndTiniestQueue& self) {
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
            [](SweepAndTiniestQueue& self) {
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
            [](SweepAndTiniestQueue& self) {
                std::vector<EdgeFaceCandidate> candidates;
                self.detect_edge_face_candidates(candidates);
                return candidates;
            },
            R"ipc_Qu8mg5v7(
            Find the candidate edge-face intersections.

            Returns:
                The candidate edge-face intersections.
            )ipc_Qu8mg5v7");

#ifdef IPC_TOOLKIT_WITH_CUDA
    py::class_<SweepAndTiniestQueueGPU, CopyMeshBroadPhase>(
        m, "SweepAndTiniestQueueGPU")
        .def(py::init())
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, double>(
                &SweepAndTiniestQueueGPU::build),
            R"ipc_Qu8mg5v7(
            Build the broad phase for static collision detection.

            Parameters:
                vertices: Vertex positions
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
                &SweepAndTiniestQueueGPU::build),
            R"ipc_Qu8mg5v7(
            Build the broad phase for continuous collision detection.

            Parameters:
                vertices_t0: Starting vertex positions
                vertices_t1: Ending vertex positions
                edges: Collision mesh edges
                faces: Collision mesh faces
                inflation_radius: Radius of inflation around all elements.
            )ipc_Qu8mg5v7",
            py::arg("vertices_t0"), py::arg("vertices_t1"), py::arg("edges"),
            py::arg("faces"), py::arg("inflation_radius") = 0)
        .def("clear", &SweepAndTiniestQueueGPU::clear, "Clear any built data.")
        .def(
            "detect_edge_vertex_candidates",
            [](SweepAndTiniestQueueGPU& self) {
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
            [](SweepAndTiniestQueueGPU& self) {
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
            [](SweepAndTiniestQueueGPU& self) {
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
            [](SweepAndTiniestQueueGPU& self) {
                std::vector<EdgeFaceCandidate> candidates;
                self.detect_edge_face_candidates(candidates);
                return candidates;
            },
            R"ipc_Qu8mg5v7(
            Find the candidate edge-face intersections.

            Returns:
                The candidate edge-face intersections.
            )ipc_Qu8mg5v7");
#endif
}
