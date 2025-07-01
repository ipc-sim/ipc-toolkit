#include <common.hpp>

#include <ipc/broad_phase/broad_phase.hpp>
#include <ipc/candidates/candidates.hpp>

using namespace ipc;

class PyBroadPhase : public BroadPhase {
public:
    using BroadPhase::BroadPhase; // Inherit constructors

    std::string name() const override
    {
        PYBIND11_OVERRIDE_PURE(std::string, BroadPhase, name);
    }

    void build(
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        const double inflation_radius = 0) override
    {
        PYBIND11_OVERRIDE(
            void, BroadPhase, build, vertices, edges, faces, inflation_radius);
    }

    void build(
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        const double inflation_radius = 0) override
    {
        PYBIND11_OVERRIDE(
            void, BroadPhase, build, vertices_t0, vertices_t1, edges, faces,
            inflation_radius);
    }

    void clear() override { PYBIND11_OVERRIDE(void, BroadPhase, clear); }

    void detect_vertex_vertex_candidates(
        std::vector<VertexVertexCandidate>& candidates) const override
    {
        py::gil_scoped_acquire gil; // Acquire GIL before calling Python code
        py::function override =
            py::get_override(this, "detect_vertex_vertex_candidates");
        if (override) {
            candidates = override().cast<std::vector<VertexVertexCandidate>>();
            return;
        }
        throw std::runtime_error(
            "Tried to call pure virtual function \"BroadPhase::detect_vertex_vertex_candidates\"");
    }

    void detect_edge_vertex_candidates(
        std::vector<EdgeVertexCandidate>& candidates) const override
    {
        py::gil_scoped_acquire gil; // Acquire GIL before calling Python code
        py::function override =
            py::get_override(this, "detect_edge_vertex_candidates");
        if (override) {
            candidates = override().cast<std::vector<EdgeVertexCandidate>>();
            return;
        }
        throw std::runtime_error(
            "Tried to call pure virtual function \"BroadPhase::detect_edge_vertex_candidates\"");
    }

    void detect_edge_edge_candidates(
        std::vector<EdgeEdgeCandidate>& candidates) const override
    {
        py::gil_scoped_acquire gil; // Acquire GIL before calling Python code
        py::function override =
            py::get_override(this, "detect_edge_edge_candidates");
        if (override) {
            candidates = override().cast<std::vector<EdgeEdgeCandidate>>();
            return;
        }
        throw std::runtime_error(
            "Tried to call pure virtual function \"BroadPhase::detect_edge_edge_candidates\"");
    }

    void detect_face_vertex_candidates(
        std::vector<FaceVertexCandidate>& candidates) const override
    {
        py::gil_scoped_acquire gil; // Acquire GIL before calling Python code
        py::function override =
            py::get_override(this, "detect_face_vertex_candidates");
        if (override) {
            candidates = override().cast<std::vector<FaceVertexCandidate>>();
            return;
        }
        throw std::runtime_error(
            "Tried to call pure virtual function \"BroadPhase::detect_face_vertex_candidates\"");
    }

    void detect_edge_face_candidates(
        std::vector<EdgeFaceCandidate>& candidates) const override
    {
        py::gil_scoped_acquire gil; // Acquire GIL before calling Python code
        py::function override =
            py::get_override(this, "detect_edge_face_candidates");
        if (override) {
            candidates = override().cast<std::vector<EdgeFaceCandidate>>();
            return;
        }
        throw std::runtime_error(
            "Tried to call pure virtual function \"BroadPhase::detect_edge_face_candidates\"");
    }

    void detect_face_face_candidates(
        std::vector<FaceFaceCandidate>& candidates) const override
    {
        py::gil_scoped_acquire gil; // Acquire GIL before calling Python code
        py::function override =
            py::get_override(this, "detect_face_face_candidates");
        if (override) {
            candidates = override().cast<std::vector<FaceFaceCandidate>>();
            return;
        }
        throw std::runtime_error(
            "Tried to call pure virtual function \"BroadPhase::detect_face_face_candidates\"");
    }
};

void define_broad_phase(py::module_& m)
{
    py::class_<BroadPhase, PyBroadPhase, std::shared_ptr<BroadPhase>>(
        m, "BroadPhase")
        .def(py::init<>())
        .def("name", &BroadPhase::name, "Get the name of the broad phase.")
        .def(
            "build",
            py::overload_cast<
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>, const double>(
                &BroadPhase::build),
            R"ipc_Qu8mg5v7(
            Build the broad phase for static collision detection.

            Parameters:
                vertices: Vertex positions
                edges: Collision mesh edges
                faces: Collision mesh faces
                inflation_radius: Radius of inflation around all elements.
            )ipc_Qu8mg5v7",
            "vertices"_a, "edges"_a, "faces"_a, "inflation_radius"_a = 0)
        .def(
            "build",
            py::overload_cast<
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>, const double>(
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
            "vertices_t0"_a, "vertices_t1"_a, "edges"_a, "faces"_a,
            "inflation_radius"_a = 0)
        .def("clear", &BroadPhase::clear, "Clear any built data.")
        .def(
            "detect_vertex_vertex_candidates",
            [](const BroadPhase& self) {
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
            [](const BroadPhase& self) {
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
            [](const BroadPhase& self) {
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
            [](const BroadPhase& self) {
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
            [](const BroadPhase& self) {
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
            "detect_face_face_candidates",
            [](const BroadPhase& self) {
                std::vector<FaceFaceCandidate> candidates;
                self.detect_face_face_candidates(candidates);
                return candidates;
            },
            R"ipc_Qu8mg5v7(
            Find the candidate face-face collisions.

            Returns:
                The candidate face-face collisions.
            )ipc_Qu8mg5v7")
        .def(
            "detect_collision_candidates",
            [](const BroadPhase& self, int dim) {
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
            "dim"_a)
        .def_readwrite(
            "can_vertices_collide", &BroadPhase::can_vertices_collide,
            "Function for determining if two vertices can collide.");
}
