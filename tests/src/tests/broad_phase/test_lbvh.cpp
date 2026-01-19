#include <ipc/config.hpp>

#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <ipc/broad_phase/bvh.hpp>
#include <ipc/broad_phase/lbvh.hpp>
// #include <ipc/broad_phase/vulkan/shaders/single_radixsort.hpp>
#include <ipc/utils/profiler.hpp>

#include <tbb/parallel_sort.h>

#include <catch2/catch_test_macros.hpp>

#include <iostream>

#include <random>

using namespace ipc;

namespace {

bool is_aabb_union(LBVH::Node parent, LBVH::Node childA, LBVH::Node childB)
{
    AABB children;
    children.min = childA.aabb_min.min(childB.aabb_min).cast<double>();
    children.max = childA.aabb_max.max(childB.aabb_max).cast<double>();
    constexpr float EPS = 1e-4f;
    return (abs(parent.aabb_max.cast<double>() - children.max) < EPS).all()
        && (abs(parent.aabb_min.cast<double>() - children.min) < EPS).all();
}

void traverse_lbvh(
    const std::vector<LBVH::Node>& lbvh_nodes,
    const uint32_t index,
    std::vector<bool>& visited)
{
    const LBVH::Node& node = lbvh_nodes[index];
    CHECK(node.is_valid());

    if (node.left == LBVH::Node::INVALID_POINTER) {
        // leaf
        CHECK(!visited[index]);
        visited[index] = true;
    } else {
        // inner node
        CHECK(!visited[index]);
        visited[index] = true;

        // verify aabbs
        LBVH::Node childA = lbvh_nodes[node.left];
        LBVH::Node childB = lbvh_nodes[node.right];

        {
            CAPTURE(
                index, node.left, node.right, node.aabb_min.transpose(),
                childA.aabb_min.transpose(), childB.aabb_min.transpose(),
                node.aabb_max.transpose(), childA.aabb_max.transpose(),
                childB.aabb_max.transpose());
            CHECK(is_aabb_union(node, childA, childB));
        }

        // continue traversal
        traverse_lbvh(lbvh_nodes, node.left, visited);
        traverse_lbvh(lbvh_nodes, node.right, visited);
    }
}

void check_valid_lbvh_nodes(const std::vector<LBVH::Node>& lbvh_nodes)
{
    std::vector<bool> visited(lbvh_nodes.size(), false);
    traverse_lbvh(lbvh_nodes, 0, visited);
    REQUIRE(
        std::all_of(visited.begin(), visited.end(), [](bool v) { return v; }));
}
} // namespace

TEST_CASE("LBVH::build", "[broad_phase][lbvh]")
{
    constexpr double inflation_radius = 1e-3;

    const std::shared_ptr<LBVH> lbvh = std::make_shared<LBVH>();

    SECTION("Static")
    {
        const std::string mesh = GENERATE("cube.ply", "bunny.ply");

        Eigen::MatrixXd vertices;
        Eigen::MatrixXi edges, faces;
        REQUIRE(tests::load_mesh(mesh, vertices, edges, faces));

        lbvh->build(vertices, edges, faces, inflation_radius);
    }
    SECTION("Dynamic")
    {
        const std::string mesh_t0 = "cloth_ball92.ply";
        const std::string mesh_t1 = "cloth_ball93.ply";

        Eigen::MatrixXd vertices_t0, vertices_t1;
        Eigen::MatrixXi edges, faces;
        REQUIRE(tests::load_mesh(mesh_t0, vertices_t0, edges, faces));
        REQUIRE(tests::load_mesh(mesh_t1, vertices_t1, edges, faces));

        lbvh->build(vertices_t0, vertices_t1, edges, faces, inflation_radius);
    }

    // -- TODO: Check the morton codes ----------------------------------------
    // -- TODO: Check the morton codes are sorted -----------------------------

    // -- Check the LBVH nodes are all reachable and contain their children ---

    check_valid_lbvh_nodes(lbvh->vertex_nodes());
    // check_valid_lbvh_nodes(lbvh->edge_nodes());
    // check_valid_lbvh_nodes(lbvh->face_nodes());

    // -- Check clear() works -------------------------------------------------
    lbvh->clear();

    // CHECK(lbvh->vertex_boxes().empty());
    // CHECK(lbvh->edge_boxes().empty());
    // CHECK(lbvh->face_boxes().empty());

    CHECK(lbvh->vertex_nodes().empty());
    CHECK(lbvh->edge_nodes().empty());
    CHECK(lbvh->face_nodes().empty());
}

namespace {

/// @brief Checks if every candidate in the expected vector is present in the actual vector.
/// @tparam Candidate The type of the candidate (e.g., VertexVertexCandidate, EdgeEdgeCandidate).
/// @param actual The vector of candidates found by the algorithm.
/// @param expected The vector of candidates that are expected to be found.
/// @return true If all candidates in `expected` are found in `actual`.
/// @return false Otherwise.
template <typename Candidate>
bool contains_all_candidates(
    std::vector<Candidate> actual, std::vector<Candidate> expected)
{
    // 1. Sort the actual candidates to prepare for set operations
    tbb::parallel_sort(actual.begin(), actual.end());
    // Ensure 'actual' has no duplicates to treat it as a mathematical set
    REQUIRE(std::unique(actual.begin(), actual.end()) == actual.end());

    // 2. Sort the expected candidates
    tbb::parallel_sort(expected.begin(), expected.end());
    // Ensure 'expected' has no duplicates to treat it as a mathematical set
    REQUIRE(std::unique(expected.begin(), expected.end()) == expected.end());

    // 3. Check if 'expected' is a subset of 'actual'
    // std::includes requires both ranges to be sorted.
    // It returns true if every element in the second range is found in the
    // first range.
    return std::includes(
        actual.begin(), actual.end(), expected.begin(), expected.end());
}

} // namespace

TEST_CASE("LBVH::detect_*_candidates", "[broad_phase][lbvh]")
{
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
    // SECTION("Puffer-Ball")
    // {
    //     mesh_t0 = "puffer-ball/20.ply";
    //     mesh_t1 = "puffer-ball/21.ply";
    // }

    Eigen::MatrixXd vertices_t0, vertices_t1;
    Eigen::MatrixXi edges, faces;
    REQUIRE(tests::load_mesh(mesh_t0, vertices_t0, edges, faces));
    REQUIRE(tests::load_mesh(mesh_t1, vertices_t1, edges, faces));

    const std::shared_ptr<LBVH> lbvh = std::make_shared<LBVH>();
    lbvh->build(vertices_t0, vertices_t1, edges, faces, inflation_radius);

    const std::shared_ptr<BVH> bvh = std::make_shared<BVH>();
    bvh->build(vertices_t0, vertices_t1, edges, faces, inflation_radius);

    // detect_vertex_vertex_candidates
    {
        std::vector<VertexVertexCandidate> vv_candidates;
        lbvh->detect_vertex_vertex_candidates(vv_candidates);

        std::vector<VertexVertexCandidate> expected_vv_candidates;
        bvh->detect_vertex_vertex_candidates(expected_vv_candidates);

        CHECK(vv_candidates.size() >= expected_vv_candidates.size());
        contains_all_candidates(vv_candidates, expected_vv_candidates);
    }

    {
        std::vector<EdgeVertexCandidate> ev_candidates;
        lbvh->detect_edge_vertex_candidates(ev_candidates);

        std::vector<EdgeVertexCandidate> expected_ev_candidates;
        bvh->detect_edge_vertex_candidates(expected_ev_candidates);

        CHECK(ev_candidates.size() >= expected_ev_candidates.size());
        contains_all_candidates(ev_candidates, expected_ev_candidates);
    }

    {
        std::vector<EdgeEdgeCandidate> ee_candidates;
        lbvh->detect_edge_edge_candidates(ee_candidates);

        std::vector<EdgeEdgeCandidate> expected_ee_candidates;
        bvh->detect_edge_edge_candidates(expected_ee_candidates);

        CHECK(ee_candidates.size() >= expected_ee_candidates.size());
        contains_all_candidates(ee_candidates, expected_ee_candidates);
    }

    {
        std::vector<FaceVertexCandidate> fv_candidates;
        lbvh->detect_face_vertex_candidates(fv_candidates);

        std::vector<FaceVertexCandidate> expected_fv_candidates;
        bvh->detect_face_vertex_candidates(expected_fv_candidates);

        CHECK(fv_candidates.size() >= expected_fv_candidates.size());
        contains_all_candidates(fv_candidates, expected_fv_candidates);
    }

    {
        std::vector<EdgeFaceCandidate> ef_candidates;
        lbvh->detect_edge_face_candidates(ef_candidates);

        std::vector<EdgeFaceCandidate> expected_ef_candidates;
        bvh->detect_edge_face_candidates(expected_ef_candidates);

        CHECK(ef_candidates.size() >= expected_ef_candidates.size());
        contains_all_candidates(ef_candidates, expected_ef_candidates);
    }

    {
        std::vector<FaceFaceCandidate> ff_candidates;
        lbvh->detect_face_face_candidates(ff_candidates);

        std::vector<FaceFaceCandidate> expected_ff_candidates;
        bvh->detect_face_face_candidates(expected_ff_candidates);

        CHECK(ff_candidates.size() >= expected_ff_candidates.size());
        contains_all_candidates(ff_candidates, expected_ff_candidates);
    }

#ifdef IPC_TOOLKIT_WITH_PROFILER
    ipc::profiler().print();
    ipc::profiler().clear();
#endif
}
