#include <ipc/config.hpp>

#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <tests/broad_phase/lbvh_validation.hpp>

#include <ipc/broad_phase/brute_force.hpp>
#include <ipc/broad_phase/spatial_hash.hpp>
#include <ipc/broad_phase/lbvh.hpp>
#include <ipc/utils/profiler.hpp>

#include <igl/edges.h>

#include <tbb/parallel_sort.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

using namespace ipc;
using ipc::tests::check_valid_lbvh_nodes;

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
    check_valid_lbvh_nodes(lbvh->edge_nodes());
    check_valid_lbvh_nodes(lbvh->face_nodes());

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

    const std::shared_ptr<SpatialHash> spatial_hash =
        std::make_shared<SpatialHash>();
    spatial_hash->build(
        vertices_t0, vertices_t1, edges, faces, inflation_radius);

    // detect_vertex_vertex_candidates
    {
        std::vector<VertexVertexCandidate> vv_candidates;
        lbvh->detect_vertex_vertex_candidates(vv_candidates);

        std::vector<VertexVertexCandidate> expected_vv_candidates;
        spatial_hash->detect_vertex_vertex_candidates(expected_vv_candidates);

        CHECK(vv_candidates.size() >= expected_vv_candidates.size());
        CHECK(contains_all_candidates(vv_candidates, expected_vv_candidates));
    }

    {
        std::vector<EdgeVertexCandidate> ev_candidates;
        lbvh->detect_edge_vertex_candidates(ev_candidates);

        std::vector<EdgeVertexCandidate> expected_ev_candidates;
        spatial_hash->detect_edge_vertex_candidates(expected_ev_candidates);

        CHECK(ev_candidates.size() >= expected_ev_candidates.size());
        CHECK(contains_all_candidates(ev_candidates, expected_ev_candidates));
    }

    {
        std::vector<EdgeEdgeCandidate> ee_candidates;
        lbvh->detect_edge_edge_candidates(ee_candidates);

        std::vector<EdgeEdgeCandidate> expected_ee_candidates;
        spatial_hash->detect_edge_edge_candidates(expected_ee_candidates);

        CHECK(ee_candidates.size() >= expected_ee_candidates.size());
        CHECK(contains_all_candidates(ee_candidates, expected_ee_candidates));
    }

    {
        std::vector<FaceVertexCandidate> fv_candidates;
        lbvh->detect_face_vertex_candidates(fv_candidates);

        std::vector<FaceVertexCandidate> expected_fv_candidates;
        spatial_hash->detect_face_vertex_candidates(expected_fv_candidates);

        CHECK(fv_candidates.size() >= expected_fv_candidates.size());
        CHECK(contains_all_candidates(fv_candidates, expected_fv_candidates));
    }

    {
        std::vector<EdgeFaceCandidate> ef_candidates;
        lbvh->detect_edge_face_candidates(ef_candidates);

        std::vector<EdgeFaceCandidate> expected_ef_candidates;
        spatial_hash->detect_edge_face_candidates(expected_ef_candidates);

        CHECK(ef_candidates.size() >= expected_ef_candidates.size());
        CHECK(contains_all_candidates(ef_candidates, expected_ef_candidates));
    }

    {
        std::vector<FaceFaceCandidate> ff_candidates;
        lbvh->detect_face_face_candidates(ff_candidates);

        std::vector<FaceFaceCandidate> expected_ff_candidates;
        spatial_hash->detect_face_face_candidates(expected_ff_candidates);

        CHECK(ff_candidates.size() >= expected_ff_candidates.size());
        CHECK(contains_all_candidates(ff_candidates, expected_ff_candidates));
    }

#ifdef IPC_TOOLKIT_WITH_PROFILER
    ipc::profiler().print();
    ipc::profiler().clear();
#endif
}

TEST_CASE("LBVH single-primitive trees", "[broad_phase][lbvh]")
{
    // A BVH over a single primitive is one node, which is both the root and a
    // leaf. When such a BVH is the traversal TARGET the descent takes a
    // dedicated branch, because the root cannot be descended into. Only two
    // detections put a BVH there that can have one node -- face-vertex (the
    // face BVH) and edge-face (the edge BVH) -- and the meshes the other tests
    // load never reduce either to a single primitive.
    //
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

    LBVH lbvh;
    lbvh.build(vertices, edges, faces, inflation_radius);

    BruteForce brute_force;
    brute_force.build(vertices, edges, faces, inflation_radius);

    // The branch under test is only reached if these really are single nodes.
    REQUIRE(lbvh.face_nodes().size() == 1);
    REQUIRE(lbvh.edge_nodes().size() == 1);

    // The LBVH rounds its AABBs outward to floats, so it may report a superset
    // of the exact (double-precision) brute-force set, never a subset.
    {
        std::vector<FaceVertexCandidate> fv_candidates, expected;
        lbvh.detect_face_vertex_candidates(fv_candidates);
        brute_force.detect_face_vertex_candidates(expected);

        // Without this the checks below would pass on an empty set, which is
        // exactly what a broken single-node branch would produce.
        REQUIRE(!expected.empty());
        CHECK(fv_candidates.size() >= expected.size());
        CHECK(contains_all_candidates(fv_candidates, expected));
    }

    {
        std::vector<EdgeFaceCandidate> ef_candidates, expected;
        lbvh.detect_edge_face_candidates(ef_candidates);
        brute_force.detect_edge_face_candidates(expected);

        REQUIRE(!expected.empty());
        CHECK(ef_candidates.size() >= expected.size());
        CHECK(contains_all_candidates(ef_candidates, expected));
    }

    // The remaining types traverse multi-node targets here, but are cheap to
    // check on a mesh this small.
    {
        std::vector<VertexVertexCandidate> vv_candidates, expected;
        lbvh.detect_vertex_vertex_candidates(vv_candidates);
        brute_force.detect_vertex_vertex_candidates(expected);
        CHECK(contains_all_candidates(vv_candidates, expected));
    }

    {
        std::vector<EdgeVertexCandidate> ev_candidates, expected;
        lbvh.detect_edge_vertex_candidates(ev_candidates);
        brute_force.detect_edge_vertex_candidates(expected);
        REQUIRE(!expected.empty());
        CHECK(contains_all_candidates(ev_candidates, expected));
    }

    // The device build reaches the same branch through the same shared code;
    // its parity with these trees is checked in test_gpu_lbvh.cu, where a
    // missing GPU skips rather than fails.
}

// A planar mesh with no inflation has zero-width vertex boxes along z, so the
// Morton normalization domain is degenerate on that axis. The reciprocal width
// must be 0 there, not infinity: 0 * inf is NaN, and converting NaN to an
// integer code is undefined. The tree must still be valid and complete.
TEST_CASE("LBVH degenerate domain", "[broad_phase][lbvh]")
{
    // A 4x4 grid of vertices in the z = 0 plane, triangulated.
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

    // Inflation 0 keeps the z extent of every box exactly zero.
    LBVH lbvh;
    lbvh.build(vertices, edges, faces, /*inflation_radius=*/0);

    check_valid_lbvh_nodes(lbvh.vertex_nodes());
    check_valid_lbvh_nodes(lbvh.edge_nodes());
    check_valid_lbvh_nodes(lbvh.face_nodes());

    // And the result still matches brute force (a superset, as usual).
    BruteForce brute_force;
    brute_force.build(vertices, edges, faces, 0);
    {
        std::vector<EdgeEdgeCandidate> candidates, expected;
        lbvh.detect_edge_edge_candidates(candidates);
        brute_force.detect_edge_edge_candidates(expected);
        REQUIRE(!expected.empty()); // coplanar neighbors do overlap
        CHECK(contains_all_candidates(candidates, expected));
    }
    {
        std::vector<FaceVertexCandidate> candidates, expected;
        lbvh.detect_face_vertex_candidates(candidates);
        brute_force.detect_face_vertex_candidates(expected);
        REQUIRE(!expected.empty());
        CHECK(contains_all_candidates(candidates, expected));
    }
}

TEST_CASE(
    "Benchmark LBVH::detect_edge_edge_candidates",
    "[!benchmark][broad_phase][lbvh]")
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
    SECTION("Puffer-Ball")
    {
        mesh_t0 = "puffer-ball/20.ply";
        mesh_t1 = "puffer-ball/21.ply";
    }

    Eigen::MatrixXd vertices_t0, vertices_t1;
    Eigen::MatrixXi edges, faces;
    REQUIRE(tests::load_mesh(mesh_t0, vertices_t0, edges, faces));
    REQUIRE(tests::load_mesh(mesh_t1, vertices_t1, edges, faces));

    const std::shared_ptr<LBVH> lbvh = std::make_shared<LBVH>();
    lbvh->build(vertices_t0, vertices_t1, edges, faces, inflation_radius);

    BENCHMARK("LBVH::detect_edge_edge_candidates")
    {
        std::vector<EdgeEdgeCandidate> ee_candidates;
        lbvh->detect_edge_edge_candidates(ee_candidates);
        return ee_candidates.size();
    };
    // The cuda::LBVH counterpart lives in test_gpu_lbvh.cu, behind the GPU
    // skip.
}

TEST_CASE("Benchmark LBVH::build", "[!benchmark][broad_phase][lbvh]")
{
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

        const std::shared_ptr<LBVH> lbvh = std::make_shared<LBVH>();

        BENCHMARK(fmt::format("LBVH::build [{}]", scene))
        {
            lbvh->build(
                vertices_t0, vertices_t1, edges, faces, inflation_radius);
            return lbvh->edge_nodes().size();
        };
        // The cuda::LBVH counterpart lives in test_gpu_lbvh.cu, behind the
        // GPU skip.
    }
}