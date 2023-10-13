#include <common.hpp>

#include <ipc/broad_phase/broad_phase.hpp>
#include <ipc/candidates/candidates.hpp>

namespace py = pybind11;
using namespace ipc;

void define_broad_phase(py::module_& m)
{
    py::enum_<BroadPhaseMethod>(
        m, "BroadPhaseMethod",
        "Enumeration of implemented broad phase methods.")
        .value("BRUTE_FORCE", BroadPhaseMethod::BRUTE_FORCE, "")
        .value("HASH_GRID", BroadPhaseMethod::HASH_GRID, "")
        .value("SPATIAL_HASH", BroadPhaseMethod::SPATIAL_HASH, "")
        .value(
            "SWEEP_AND_TINIEST_QUEUE",
            BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE, "")
        .value(
            "SWEEP_AND_TINIEST_QUEUE_GPU",
            BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE_GPU, "")
        .export_values();

    py::class_<BroadPhase>(m, "BroadPhase")
        .def_static(
            "make_broad_phase", &BroadPhase::make_broad_phase,
            R"ipc_Qu8mg5v7(
            Construct a registered broad phase object.

            Parameters:
                broad_phase_method: The broad phase method to use.

            Returns:
                The constructed broad phase object.
            )ipc_Qu8mg5v7",
            py::arg("broad_phase_method"))
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, double>(&BroadPhase::build),
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
                &BroadPhase::build),
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
        .def("clear", &BroadPhase::clear, "Clear any built data.")
        .def(
            "detect_edge_vertex_candidates",
            [](BroadPhase& self) {
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
            [](BroadPhase& self) {
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
            [](BroadPhase& self) {
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
            [](BroadPhase& self) {
                std::vector<EdgeFaceCandidate> candidates;
                self.detect_edge_face_candidates(candidates);
                return candidates;
            },
            R"ipc_Qu8mg5v7(
            Find the candidate edge-face intersections.

            Returns:
                The candidate edge-face intersections.
            )ipc_Qu8mg5v7")
        .def(
            "detect_collision_candidates",
            [](BroadPhase& self, int dim) {
                Candidates candidates;
                self.detect_collision_candidates(dim, candidates);
                return candidates;
            },
            R"ipc_Qu8mg5v7(
            Detect all collision candidates needed for a given dimensional simulation.

            Parameters:
                dim: The dimension of the simulation (i.e., 2 or 3).
            )ipc_Qu8mg5v7",
            py::arg("dim"))
        .def_readwrite(
            "can_vertices_collide", &BroadPhase::can_vertices_collide,
            "Function for determining if two vertices can collide.");
}
