#include <common.hpp>

#include <ipc/broad_phase/brute_force.hpp>

namespace py = pybind11;
using namespace ipc;

void define_brute_force(py::module_& m)
{
    py::class_<BruteForce, BroadPhase>(m, "BruteForce")
        .def(
            "detect_vertex_vertex_candidates",
            [](BruteForce& self) {
                std::vector<VertexVertexCandidate> candidates;
                self.detect_vertex_vertex_candidates(candidates);
                return candidates;
            },
            "Find the candidate vertex-vertex collisions.")
        .def(
            "detect_edge_vertex_candidates",
            [](BruteForce& self) {
                std::vector<EdgeVertexCandidate> candidates;
                self.detect_edge_vertex_candidates(candidates);
                return candidates;
            },
            "Find the candidate edge-vertex collisions.")
        .def(
            "detect_edge_edge_candidates",
            [](BruteForce& self) {
                std::vector<EdgeEdgeCandidate> candidates;
                self.detect_edge_edge_candidates(candidates);
                return candidates;
            },
            "Find the candidate edge-edge collisions.")
        .def(
            "detect_face_vertex_candidates",
            [](BruteForce& self) {
                std::vector<FaceVertexCandidate> candidates;
                self.detect_face_vertex_candidates(candidates);
                return candidates;
            },
            "Find the candidate face-vertex collisions.")
        .def(
            "detect_edge_face_candidates",
            [](BruteForce& self) {
                std::vector<EdgeFaceCandidate> candidates;
                self.detect_edge_face_candidates(candidates);
                return candidates;
            },
            "Find the candidate edge-face intersections.");
}
