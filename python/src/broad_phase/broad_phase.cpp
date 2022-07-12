#include "../common.hpp"

#include <ipc/broad_phase/broad_phase.hpp>

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
#ifdef IPC_TOOLKIT_WITH_CUDA
        .value(
            "SWEEP_AND_TINIEST_QUEUE_GPU",
            BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE_GPU, "")
#endif
        .export_values();

    py::class_<BroadPhase>(m, "BroadPhase")
        .def_static(
            "make_broad_phase", &BroadPhase::make_broad_phase,
            R"ipc_Qu8mg5v7(
            Construct a registered broad phase object.

            Parameters:
                method: The broad phase method to use.

            Returns:
                The constructed broad phase object.
            )ipc_Qu8mg5v7",
            py::arg("method"))
        .def(
            "build",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, double>(&BroadPhase::build),
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
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, double>(
                &BroadPhase::build),
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
                candidates: The detected collision candidates.
            )ipc_Qu8mg5v7",
            py::arg("dim"))
        .def_readwrite(
            "can_vertices_collide", &BroadPhase::can_vertices_collide,
            "Function for determining if two vertices can collide.");

    m.def(
        "construct_collision_candidates",
        [](const CollisionMesh& mesh, const Eigen::MatrixXd& V,
           double inflation_radius = 0,
           const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID) {
            Candidates candidates;
            construct_collision_candidates(
                mesh, V, candidates, inflation_radius, method);
            return candidates;
        },
        R"ipc_Qu8mg5v7(
        Construct a set of discrete collision detection candidates.

        Parameters:
            mesh: The surface of the contact mesh.
            V: Surface Vertex positions at start as rows of a matrix.
            inflation_radius: Amount to inflate the bounding boxes.
            method: Broad phase method to use.

        Returns:
            The constructed candidate set as output.
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("V"), py::arg("inflation_radius") = 0,
        py::arg("method") = BroadPhaseMethod::HASH_GRID);

    m.def(
        "construct_collision_candidates",
        [](const CollisionMesh& mesh, const Eigen::MatrixXd& V0,
           const Eigen::MatrixXd& V1, double inflation_radius = 0,
           const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID) {
            Candidates candidates;
            construct_collision_candidates(
                mesh, V0, V1, candidates, inflation_radius, method);
            return candidates;
        },
        R"ipc_Qu8mg5v7(
        Construct a set of continuous collision detection candidates.

        Note:
            Assumes the trajectory is linear.

        Parameters:
            mesh: The surface of the contact mesh.
            V0: Surface vertex positions at start as rows of a matrix.
            V1: Surface vertex positions at end as rows of a matrix.
            inflation_radius: Amount to inflate the bounding boxes.
            method: Broad phase method to use.

        Returns:
            The constructed candidate set as output.
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("V0"), py::arg("V1"),
        py::arg("inflation_radius") = 0,
        py::arg("method") = BroadPhaseMethod::HASH_GRID);
}
