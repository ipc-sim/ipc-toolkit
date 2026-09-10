// Validates the GPU-built LBVH (ipc::cuda::LBVH) against the CPU ipc::LBVH:
// ipc::cuda::LBVH builds the vertex/edge/face AABBs and BVHs entirely on the
// device. The copied-back trees must be structurally valid (every node
// reachable exactly once, every internal AABB the union of its children, leaf
// set = {0..n-1}) and must agree with the CPU build on node count and root
// AABB (an order-independent union of identically-inflated boxes). The
// detected candidate sets must be exactly equal.

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <tests/config.hpp>
#include <tests/gpu_utils.hpp>
#include <tests/utils.hpp>
#include <tests/broad_phase/lbvh_validation.hpp>

#include <ipc/broad_phase/lbvh.hpp>
#include <ipc/broad_phase/cuda/lbvh.hpp>

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <igl/edges.h>
#include <igl/readCSV.h>

#include <algorithm>
#include <array>
#include <string>
#include <vector>

using namespace ipc;

namespace {

// The GPU and CPU candidate sets are determined by the (bit-identical) box
// overlaps + the same can_*_collide predicate, independent of tree structure,
// so they must be exactly equal as sets.
template <typename Candidate>
void compare_candidates_exact(
    std::vector<Candidate> gpu, std::vector<Candidate> cpu)
{
    std::sort(gpu.begin(), gpu.end());
    std::sort(cpu.begin(), cpu.end());
    CHECK(gpu.size() == cpu.size());
    CHECK(gpu == cpu);
}

} // namespace

TEST_CASE("GPU LBVH build", "[broad_phase][lbvh][cuda][gpu]")
{
    tests::skip_if_no_cuda_device();

    constexpr double inflation_radius = 1e-3;

    const std::string mesh = GENERATE("cube.ply", "bunny.ply");
    CAPTURE(mesh);

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    REQUIRE(tests::load_mesh(mesh, vertices, edges, faces));

    // GPU build (boxes + BVHs all on the device).
    cuda::LBVH gpu_lbvh;
    gpu_lbvh.build(vertices, edges, faces, inflation_radius);

    // CPU reference.
    LBVH cpu_lbvh;
    cpu_lbvh.build(vertices, edges, faces, inflation_radius);

    LBVH::Nodes nodes;
    LBVH::RightmostLeaves rightmost;

    SECTION("vertices")
    {
        gpu_lbvh.vertex_nodes_to_host(nodes, rightmost);
        tests::check_lbvh_nodes_match(nodes, cpu_lbvh.vertex_nodes());
    }
    SECTION("edges")
    {
        gpu_lbvh.edge_nodes_to_host(nodes, rightmost);
        tests::check_lbvh_nodes_match(nodes, cpu_lbvh.edge_nodes());
    }
    SECTION("faces")
    {
        gpu_lbvh.face_nodes_to_host(nodes, rightmost);
        tests::check_lbvh_nodes_match(nodes, cpu_lbvh.face_nodes());
    }

    // clear() empties the device trees.
    gpu_lbvh.clear();
    CHECK(gpu_lbvh.num_vertex_nodes() == 0);
    gpu_lbvh.vertex_nodes_to_host(nodes, rightmost);
    CHECK(nodes.empty());
}

TEST_CASE("GPU LBVH detect candidates", "[broad_phase][lbvh][cuda][gpu]")
{
    tests::skip_if_no_cuda_device();

    constexpr double inflation_radius = 0;

    std::string mesh_t0, mesh_t1;
    SECTION("Two cubes")
    {
        mesh_t0 = "two-cubes-far.ply";
        mesh_t1 = "two-cubes-intersecting.ply";
    }
    SECTION("Cloth-Ball")
    {
        mesh_t0 = "cloth_ball92.ply";
        mesh_t1 = "cloth_ball93.ply";
    }

    Eigen::MatrixXd vertices_t0, vertices_t1;
    Eigen::MatrixXi edges, faces;
    REQUIRE(tests::load_mesh(mesh_t0, vertices_t0, edges, faces));
    REQUIRE(tests::load_mesh(mesh_t1, vertices_t1, edges, faces));

    cuda::LBVH gpu_lbvh;
    gpu_lbvh.build(vertices_t0, vertices_t1, edges, faces, inflation_radius);

    LBVH cpu_lbvh;
    cpu_lbvh.build(vertices_t0, vertices_t1, edges, faces, inflation_radius);

    {
        std::vector<VertexVertexCandidate> gpu_c, cpu_c;
        gpu_lbvh.detect_vertex_vertex_candidates(gpu_c);
        cpu_lbvh.detect_vertex_vertex_candidates(cpu_c);
        compare_candidates_exact(gpu_c, cpu_c);
    }
    {
        std::vector<EdgeVertexCandidate> gpu_c, cpu_c;
        gpu_lbvh.detect_edge_vertex_candidates(gpu_c);
        cpu_lbvh.detect_edge_vertex_candidates(cpu_c);
        compare_candidates_exact(gpu_c, cpu_c);
    }
    {
        std::vector<EdgeEdgeCandidate> gpu_c, cpu_c;
        gpu_lbvh.detect_edge_edge_candidates(gpu_c);
        cpu_lbvh.detect_edge_edge_candidates(cpu_c);
        compare_candidates_exact(gpu_c, cpu_c);

        // With the default (accept-all) filter the device-resident buffer is
        // already the exact set (no host trimming needed), and it holds the
        // same pairs the host variant materialized.
        const cuda::LBVH::DeviceCandidateView view =
            gpu_lbvh.detect_edge_edge_candidates_device();
        REQUIRE(view.size == cpu_c.size());
        std::vector<int32_t> a(view.size), b(view.size);
        REQUIRE_CUDA(cudaMemcpy(
            a.data(), view.a, view.size * sizeof(int32_t),
            cudaMemcpyDeviceToHost));
        REQUIRE_CUDA(cudaMemcpy(
            b.data(), view.b, view.size * sizeof(int32_t),
            cudaMemcpyDeviceToHost));
        std::vector<EdgeEdgeCandidate> view_c;
        for (size_t k = 0; k < view.size; ++k) {
            view_c.emplace_back(a[k], b[k]);
        }
        compare_candidates_exact(view_c, cpu_c);

        // Like every BroadPhase, detection clears its output first: a second
        // call replaces the vector rather than doubling it.
        gpu_lbvh.detect_edge_edge_candidates(gpu_c);
        CHECK(gpu_c.size() == cpu_c.size());
    }
    {
        std::vector<FaceVertexCandidate> gpu_c, cpu_c;
        gpu_lbvh.detect_face_vertex_candidates(gpu_c);
        cpu_lbvh.detect_face_vertex_candidates(cpu_c);
        compare_candidates_exact(gpu_c, cpu_c);
    }
    {
        std::vector<EdgeFaceCandidate> gpu_c, cpu_c;
        gpu_lbvh.detect_edge_face_candidates(gpu_c);
        cpu_lbvh.detect_edge_face_candidates(cpu_c);
        compare_candidates_exact(gpu_c, cpu_c);
    }
    {
        std::vector<FaceFaceCandidate> gpu_c, cpu_c;
        gpu_lbvh.detect_face_face_candidates(gpu_c);
        cpu_lbvh.detect_face_face_candidates(cpu_c);
        compare_candidates_exact(gpu_c, cpu_c);
    }
}

// Exercises the host-fallback path: a non-trivial user vertex filter is not
// device-representable yet, so the device emits the connectivity-filtered
// superset and the host trims it. The result must still match the CPU exactly.
TEST_CASE(
    "GPU LBVH detect candidates (custom filter)",
    "[broad_phase][lbvh][cuda][gpu]")
{
    tests::skip_if_no_cuda_device();

    Eigen::MatrixXd vertices_t0, vertices_t1;
    Eigen::MatrixXi edges, faces;
    REQUIRE(tests::load_mesh("two-cubes-far.ply", vertices_t0, edges, faces));
    REQUIRE(
        tests::load_mesh(
            "two-cubes-intersecting.ply", vertices_t1, edges, faces));

    // An arbitrary (not device-representable) filter -> host fallback.
    const auto filter = [](size_t a, size_t b) { return ((a + b) % 2) == 0; };

    cuda::LBVH gpu_lbvh;
    gpu_lbvh.can_vertices_collide = filter;
    gpu_lbvh.build(vertices_t0, vertices_t1, edges, faces, 0);
    REQUIRE_FALSE(gpu_lbvh.can_vertices_collide.accepts_all());

    LBVH cpu_lbvh;
    cpu_lbvh.can_vertices_collide = filter;
    cpu_lbvh.build(vertices_t0, vertices_t1, edges, faces, 0);

    {
        std::vector<EdgeEdgeCandidate> gpu_c, cpu_c;
        gpu_lbvh.detect_edge_edge_candidates(gpu_c);
        cpu_lbvh.detect_edge_edge_candidates(cpu_c);
        compare_candidates_exact(gpu_c, cpu_c);

        // The device view is the connectivity-filtered superset the host
        // trimmed: never smaller than the exact set.
        CHECK(
            gpu_lbvh.detect_edge_edge_candidates_device().size >= cpu_c.size());
    }
    {
        std::vector<FaceVertexCandidate> gpu_c, cpu_c;
        gpu_lbvh.detect_face_vertex_candidates(gpu_c);
        cpu_lbvh.detect_face_vertex_candidates(cpu_c);
        compare_candidates_exact(gpu_c, cpu_c);
    }
    {
        std::vector<EdgeFaceCandidate> gpu_c, cpu_c;
        gpu_lbvh.detect_edge_face_candidates(gpu_c);
        cpu_lbvh.detect_edge_face_candidates(cpu_c);
        compare_candidates_exact(gpu_c, cpu_c);
    }
}

// 2D input has no faces; ipc::AABB zero-pads the unused z component without
// inflating it (see build_vertex_boxes_kernel in lbvh.cu), so this also
// exercises that padding path against the CPU's exact behavior.
TEST_CASE("GPU LBVH 2D build and detect", "[broad_phase][lbvh][cuda][gpu]")
{
    tests::skip_if_no_cuda_device();

    Eigen::MatrixXd tmp;
    REQUIRE(igl::readCSV((tests::DATA_DIR / "mesh-2D/V_t0.csv").string(), tmp));
    const Eigen::MatrixXd V0 = tmp.leftCols(2);
    REQUIRE(igl::readCSV((tests::DATA_DIR / "mesh-2D/V_t1.csv").string(), tmp));
    const Eigen::MatrixXd V1 = tmp.leftCols(2);
    Eigen::MatrixXi E;
    REQUIRE(igl::readCSV((tests::DATA_DIR / "mesh-2D/E.csv").string(), E));
    E.array() -= 1; // Convert from OBJ format to 0-indexed
    const Eigen::MatrixXi F(0, 3);

    constexpr double inflation_radius = 1e-3;

    cuda::LBVH gpu_lbvh;
    gpu_lbvh.build(V0, V1, E, F, inflation_radius);

    LBVH cpu_lbvh;
    cpu_lbvh.build(V0, V1, E, F, inflation_radius);

    // -- Build parity (structure + root AABB, same checks as the 3D case). --
    LBVH::Nodes nodes;
    LBVH::RightmostLeaves rightmost;
    gpu_lbvh.vertex_nodes_to_host(nodes, rightmost);
    tests::check_lbvh_nodes_match(nodes, cpu_lbvh.vertex_nodes());
    gpu_lbvh.edge_nodes_to_host(nodes, rightmost);
    tests::check_lbvh_nodes_match(nodes, cpu_lbvh.edge_nodes());
    CHECK(gpu_lbvh.num_face_nodes() == 0);

    // -- Detection parity (only edge-vertex is meaningful in 2D; mirrors
    //    BroadPhase::detect_collision_candidates's dim == 2 branch). --
    std::vector<EdgeVertexCandidate> gpu_c, cpu_c;
    gpu_lbvh.detect_edge_vertex_candidates(gpu_c);
    cpu_lbvh.detect_edge_vertex_candidates(cpu_c);
    compare_candidates_exact(gpu_c, cpu_c);
    CHECK(!gpu_c.empty());
}

// A BVH over a single primitive is one node, both root and leaf, and the
// shared descent takes a dedicated branch for such a TARGET. The CPU half of
// this scenario is checked against brute force in test_lbvh.cpp; here the
// device build -- which reaches the same branch through the same shared code
// but from a kernel -- must agree with the host on exactly these trees.
TEST_CASE("GPU LBVH single-primitive trees", "[broad_phase][lbvh][cuda][gpu]")
{
    tests::skip_if_no_cuda_device();

    // One face and one edge, sharing no vertices so the connectivity filter
    // keeps the pair, and inflated enough that the AABBs actually overlap.
    Eigen::MatrixXd vertices(5, 3);
    vertices << 0.00, 0.00, 0.00, // 0 |
        1.00, 0.00, 0.00,         // 1 |- the face
        0.00, 1.00, 0.00,         // 2 |
        0.05, 0.05, 0.05,         // 3 |- the edge
        0.15, 0.05, 0.05;         // 4 |

    Eigen::MatrixXi edges(1, 2);
    edges << 3, 4;

    Eigen::MatrixXi faces(1, 3);
    faces << 0, 1, 2;

    constexpr double inflation_radius = 0.1;

    cuda::LBVH gpu_lbvh;
    gpu_lbvh.build(vertices, edges, faces, inflation_radius);

    LBVH cpu_lbvh;
    cpu_lbvh.build(vertices, edges, faces, inflation_radius);

    // The branch under test is only reached if these really are single nodes.
    REQUIRE(gpu_lbvh.num_face_nodes() == 1);
    REQUIRE(gpu_lbvh.num_edge_nodes() == 1);
    REQUIRE(cpu_lbvh.face_nodes().size() == 1);
    REQUIRE(cpu_lbvh.edge_nodes().size() == 1);

    LBVH::Nodes nodes;
    LBVH::RightmostLeaves rightmost;
    gpu_lbvh.face_nodes_to_host(nodes, rightmost);
    tests::check_lbvh_nodes_match(nodes, cpu_lbvh.face_nodes());
    gpu_lbvh.edge_nodes_to_host(nodes, rightmost);
    tests::check_lbvh_nodes_match(nodes, cpu_lbvh.edge_nodes());

    {
        std::vector<FaceVertexCandidate> gpu_c, cpu_c;
        gpu_lbvh.detect_face_vertex_candidates(gpu_c);
        cpu_lbvh.detect_face_vertex_candidates(cpu_c);
        // Without this the checks would pass on an empty set, which is
        // exactly what a broken single-node branch would produce.
        REQUIRE(!cpu_c.empty());
        compare_candidates_exact(gpu_c, cpu_c);
    }
    {
        std::vector<EdgeFaceCandidate> gpu_c, cpu_c;
        gpu_lbvh.detect_edge_face_candidates(gpu_c);
        cpu_lbvh.detect_edge_face_candidates(cpu_c);
        REQUIRE(!cpu_c.empty());
        compare_candidates_exact(gpu_c, cpu_c);
    }
    {
        std::vector<EdgeVertexCandidate> gpu_c, cpu_c;
        gpu_lbvh.detect_edge_vertex_candidates(gpu_c);
        cpu_lbvh.detect_edge_vertex_candidates(cpu_c);
        REQUIRE(!cpu_c.empty());
        compare_candidates_exact(gpu_c, cpu_c);
    }
}

// A planar mesh with no inflation makes the Morton normalization domain
// zero-width along z (see morton_domain_width_inv()). The device build must
// take the same guarded path as the host and produce the same trees.
TEST_CASE("GPU LBVH degenerate domain", "[broad_phase][lbvh][cuda][gpu]")
{
    tests::skip_if_no_cuda_device();

    constexpr int N = 4;
    Eigen::MatrixXd vertices(N * N, 3);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            vertices.row(N * i + j) << i, j, 0.0;
        }
    }
    Eigen::MatrixXi faces(2 * (N - 1) * (N - 1), 3);
    for (int i = 0, f = 0; i < N - 1; ++i) {
        for (int j = 0; j < N - 1; ++j) {
            const int v00 = N * i + j, v10 = v00 + N;
            faces.row(f++) << v00, v10, v00 + 1;
            faces.row(f++) << v10, v10 + 1, v00 + 1;
        }
    }
    Eigen::MatrixXi edges;
    igl::edges(faces, edges);

    cuda::LBVH gpu_lbvh;
    gpu_lbvh.build(vertices, edges, faces, /*inflation_radius=*/0);

    LBVH cpu_lbvh;
    cpu_lbvh.build(vertices, edges, faces, 0);

    LBVH::Nodes nodes;
    LBVH::RightmostLeaves rightmost;
    gpu_lbvh.vertex_nodes_to_host(nodes, rightmost);
    tests::check_lbvh_nodes_match(nodes, cpu_lbvh.vertex_nodes());
    gpu_lbvh.edge_nodes_to_host(nodes, rightmost);
    tests::check_lbvh_nodes_match(nodes, cpu_lbvh.edge_nodes());
    gpu_lbvh.face_nodes_to_host(nodes, rightmost);
    tests::check_lbvh_nodes_match(nodes, cpu_lbvh.face_nodes());

    std::vector<EdgeEdgeCandidate> gpu_c, cpu_c;
    gpu_lbvh.detect_edge_edge_candidates(gpu_c);
    cpu_lbvh.detect_edge_edge_candidates(cpu_c);
    REQUIRE(!cpu_c.empty());
    compare_candidates_exact(gpu_c, cpu_c);
}

// The moved-from object must stay usable: it is cleared, not left with a null
// implementation.
TEST_CASE("GPU LBVH move", "[broad_phase][lbvh][cuda][gpu]")
{
    tests::skip_if_no_cuda_device();

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    REQUIRE(tests::load_mesh("cube.ply", vertices, edges, faces));

    cuda::LBVH a;
    a.build(vertices, edges, faces, 1e-3);
    const size_t num_nodes = a.num_face_nodes();
    REQUIRE(num_nodes > 1);

    cuda::LBVH b(std::move(a));
    CHECK(b.num_face_nodes() == num_nodes);
    CHECK(a.num_face_nodes() == 0); // NOLINT(bugprone-use-after-move)

    std::vector<FaceFaceCandidate> candidates;
    a.detect_face_face_candidates(candidates); // cleared, not null
    CHECK(candidates.empty());

    a.build(vertices, edges, faces, 1e-3); // and rebuildable
    CHECK(a.num_face_nodes() == num_nodes);

    cuda::LBVH c;
    c = std::move(b);
    CHECK(c.num_face_nodes() == num_nodes);
    CHECK(b.num_face_nodes() == 0); // NOLINT(bugprone-use-after-move)
}

// ---------------------------------------------------------------------------
// Benchmarks. Hidden ([!benchmark]) and GPU-gated like every other case here,
// so a CUDA build without a device skips rather than fails them.

TEST_CASE(
    "Benchmark cuda::LBVH::detect_edge_edge_candidates",
    "[!benchmark][broad_phase][lbvh][cuda][gpu]")
{
    tests::skip_if_no_cuda_device();

    constexpr double inflation_radius = 0;

    std::string mesh_t0, mesh_t1;
    SECTION("Two cubes")
    {
        mesh_t0 = "two-cubes-far.ply";
        mesh_t1 = "two-cubes-intersecting.ply";
    }
    SECTION("Cloth-Ball")
    {
        mesh_t0 = "cloth_ball92.ply";
        mesh_t1 = "cloth_ball93.ply";
    }
#ifdef NDEBUG
    SECTION("Armadillo-Rollers")
    {
        mesh_t0 = "armadillo-rollers/326.ply";
        mesh_t1 = "armadillo-rollers/327.ply";
    }
    SECTION("Cloth-Funnel")
    {
        mesh_t0 = "cloth-funnel/227.ply";
        mesh_t1 = "cloth-funnel/228.ply";
    }
    SECTION("N-Body-Simulation")
    {
        mesh_t0 = "n-body-simulation/balls16_18.ply";
        mesh_t1 = "n-body-simulation/balls16_19.ply";
    }
    SECTION("Rod-Twist")
    {
        mesh_t0 = "rod-twist/3036.ply";
        mesh_t1 = "rod-twist/3037.ply";
    }
#endif
    SECTION("Puffer-Ball")
    {
        mesh_t0 = "puffer-ball/20.ply";
        mesh_t1 = "puffer-ball/21.ply";
    }

    Eigen::MatrixXd vertices_t0, vertices_t1;
    Eigen::MatrixXi edges, faces;
    REQUIRE(tests::load_mesh(mesh_t0, vertices_t0, edges, faces));
    REQUIRE(tests::load_mesh(mesh_t1, vertices_t1, edges, faces));

    cuda::LBVH gpu_lbvh;
    gpu_lbvh.build(vertices_t0, vertices_t1, edges, faces, inflation_radius);
    // Warm up the CUDA context so the first sample is not skewed by lazy
    // context/allocation initialization.
    {
        std::vector<EdgeEdgeCandidate> warmup;
        gpu_lbvh.detect_edge_edge_candidates(warmup);
    }

    BENCHMARK("cuda::LBVH::detect_edge_edge_candidates")
    {
        std::vector<EdgeEdgeCandidate> ee_candidates;
        gpu_lbvh.detect_edge_edge_candidates(ee_candidates);
        return ee_candidates.size();
    };
}

TEST_CASE(
    "Benchmark cuda::LBVH::build", "[!benchmark][broad_phase][lbvh][cuda][gpu]")
{
    tests::skip_if_no_cuda_device();

    constexpr double inflation_radius = 0;

    struct Scene {
        std::string name, mesh_t0, mesh_t1;
    };

#ifdef NDEBUG
    constexpr int NUM_SCENES = 6;
#else
    constexpr int NUM_SCENES = 1;
#endif

    const std::array<Scene, NUM_SCENES> scenes = { {
        Scene { "Cloth-Ball", "cloth_ball92.ply", "cloth_ball93.ply" },
#ifdef NDEBUG
        Scene { "Cloth-Funnel", "cloth-funnel/227.ply",
                "cloth-funnel/228.ply" },
        Scene { "Armadillo-Rollers", "armadillo-rollers/326.ply",
                "armadillo-rollers/327.ply" },
        Scene { "Rod-Twist", "rod-twist/3036.ply", "rod-twist/3037.ply" },
        Scene { "N-Body-Simulation", "n-body-simulation/balls16_18.ply",
                "n-body-simulation/balls16_19.ply" },
        Scene { "Puffer-Ball", "puffer-ball/20.ply", "puffer-ball/21.ply" },
#endif
    } };

    for (const auto& [scene, mesh_t0, mesh_t1] : scenes) {
        Eigen::MatrixXd vertices_t0, vertices_t1;
        Eigen::MatrixXi edges, faces;
        REQUIRE(tests::load_mesh(mesh_t0, vertices_t0, edges, faces));
        REQUIRE(tests::load_mesh(mesh_t1, vertices_t1, edges, faces));

        cuda::LBVH gpu_lbvh;
        // Warm up the CUDA context so the first sample is not skewed by lazy
        // context/allocation initialization.
        gpu_lbvh.build(
            vertices_t0, vertices_t1, edges, faces, inflation_radius);

        BENCHMARK("cuda::LBVH::build [" + scene + "]")
        {
            gpu_lbvh.build(
                vertices_t0, vertices_t1, edges, faces, inflation_radius);
            return gpu_lbvh.num_edge_nodes();
        };
    }
}

#endif // IPC_TOOLKIT_WITH_CUDA
