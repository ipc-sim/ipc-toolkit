// Validates the GPU-built LBVH (ipc::cuda::LBVH) against the CPU ipc::LBVH:
// ipc::cuda::LBVH builds the vertex/edge/face AABBs and BVHs entirely on the
// device. The copied-back trees must be structurally valid (every node
// reachable exactly once, every internal AABB the union of its children, leaf
// set = {0..n-1}) and must agree with the CPU build on node count and root
// AABB (an order-independent union of identically-inflated boxes).

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <ipc/broad_phase/lbvh.hpp>
#include <ipc/broad_phase/cuda/lbvh.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <cuda_runtime.h>
#include <igl/readCSV.h>

#include <algorithm>
#include <vector>

using namespace ipc;

namespace {

bool has_cuda_device()
{
    int n = 0;
    return cudaGetDeviceCount(&n) == cudaSuccess && n > 0;
}

bool is_aabb_union(
    const LBVH::Node& parent,
    const LBVH::Node& child_a,
    const LBVH::Node& child_b)
{
    const Eigen::Array3d cmin =
        child_a.aabb_min.min(child_b.aabb_min).cast<double>();
    const Eigen::Array3d cmax =
        child_a.aabb_max.max(child_b.aabb_max).cast<double>();
    constexpr float EPS = 1e-4f;
    return (abs(parent.aabb_max.cast<double>() - cmax) < EPS).all()
        && (abs(parent.aabb_min.cast<double>() - cmin) < EPS).all();
}

// Recursively verify reachability (each node visited exactly once) and that
// every internal node's AABB is the union of its children's. Collects the leaf
// primitive ids that are reached.
void traverse_and_check(
    const LBVH::Nodes& nodes,
    const int32_t index,
    std::vector<bool>& visited,
    std::vector<int32_t>& reached_leaves)
{
    REQUIRE(index >= 0);
    REQUIRE(index < int32_t(nodes.size()));
    const LBVH::Node& node = nodes[index];
    CHECK(node.is_valid());
    CHECK(!visited[index]);
    visited[index] = true;

    if (node.is_leaf()) {
        reached_leaves.push_back(node.primitive_id);
        return;
    }

    const LBVH::Node& child_a = nodes[node.left];
    const LBVH::Node& child_b = nodes[node.right];
    {
        CAPTURE(index, node.left, node.right);
        CHECK(is_aabb_union(node, child_a, child_b));
    }
    traverse_and_check(nodes, node.left, visited, reached_leaves);
    traverse_and_check(nodes, node.right, visited, reached_leaves);
}

// Validate one device-built tree (copied back to the host) against the
// corresponding CPU-built node array.
void check_tree(const LBVH::Nodes& nodes, const LBVH::Nodes& cpu_nodes)
{
    if (nodes.size() <= 1) {
        return; // single-node trees are not exercised here
    }
    REQUIRE(nodes.size() == cpu_nodes.size());
    REQUIRE(nodes.size() % 2 == 1); // 2n - 1
    const size_t n_leaves = (nodes.size() + 1) / 2;

    // -- Structural validity: reachable-once + AABB unions. --
    std::vector<bool> visited(nodes.size(), false);
    std::vector<int32_t> reached_leaves;
    traverse_and_check(nodes, 0, visited, reached_leaves);
    CHECK(
        std::all_of(visited.begin(), visited.end(), [](bool v) { return v; }));

    // -- Leaf set must be exactly {0, ..., n_leaves - 1}. --
    REQUIRE(reached_leaves.size() == n_leaves);
    std::sort(reached_leaves.begin(), reached_leaves.end());
    for (size_t i = 0; i < reached_leaves.size(); ++i) {
        CHECK(reached_leaves[i] == int32_t(i));
    }

    // -- Root AABB must equal the CPU root AABB (an order-independent union of
    //    identically-inflated boxes). --
    constexpr float EPS = 1e-4f;
    CHECK((abs(nodes[0].aabb_min.cast<double>()
               - cpu_nodes[0].aabb_min.cast<double>())
           < EPS)
              .all());
    CHECK((abs(nodes[0].aabb_max.cast<double>()
               - cpu_nodes[0].aabb_max.cast<double>())
           < EPS)
              .all());
}

} // namespace

TEST_CASE("GPU LBVH build", "[broad_phase][lbvh][cuda][gpu]")
{
    if (!has_cuda_device()) {
        SKIP("No CUDA device available; kernels compiled but not executed.");
    }

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
        check_tree(nodes, cpu_lbvh.vertex_nodes());
    }
    SECTION("edges")
    {
        gpu_lbvh.edge_nodes_to_host(nodes, rightmost);
        check_tree(nodes, cpu_lbvh.edge_nodes());
    }
    SECTION("faces")
    {
        gpu_lbvh.face_nodes_to_host(nodes, rightmost);
        check_tree(nodes, cpu_lbvh.face_nodes());
    }

    // clear() empties the device trees.
    gpu_lbvh.clear();
    gpu_lbvh.vertex_nodes_to_host(nodes, rightmost);
    CHECK(nodes.empty());
}

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

TEST_CASE("GPU LBVH detect candidates", "[broad_phase][lbvh][cuda][gpu]")
{
    if (!has_cuda_device()) {
        SKIP("No CUDA device available; kernels compiled but not executed.");
    }

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
        // already the exact set (no host trimming needed).
        CHECK(
            gpu_lbvh.detect_edge_edge_candidates_device().size == cpu_c.size());
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
    if (!has_cuda_device()) {
        SKIP("No CUDA device available; kernels compiled but not executed.");
    }

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
// inflating it (see build_vertex_boxes_{static,dynamic}_kernel in lbvh.cu), so
// this also exercises that padding path against the CPU's exact behavior.
TEST_CASE("GPU LBVH 2D build and detect", "[broad_phase][lbvh][cuda][gpu]")
{
    if (!has_cuda_device()) {
        SKIP("No CUDA device available; kernels compiled but not executed.");
    }

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
    check_tree(nodes, cpu_lbvh.vertex_nodes());
    gpu_lbvh.edge_nodes_to_host(nodes, rightmost);
    check_tree(nodes, cpu_lbvh.edge_nodes());

    // -- Detection parity (only edge-vertex is meaningful in 2D; mirrors
    //    BroadPhase::detect_collision_candidates's dim == 2 branch). --
    std::vector<EdgeVertexCandidate> gpu_c, cpu_c;
    gpu_lbvh.detect_edge_vertex_candidates(gpu_c);
    cpu_lbvh.detect_edge_vertex_candidates(cpu_c);
    compare_candidates_exact(gpu_c, cpu_c);
    CHECK(!gpu_c.empty());
}

#endif // IPC_TOOLKIT_WITH_CUDA
