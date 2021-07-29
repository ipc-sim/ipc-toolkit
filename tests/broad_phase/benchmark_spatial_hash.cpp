#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <igl/IO>
#include <igl/edges.h>
#include <igl/Timer.h>

#include <ipc/ipc.hpp>
#include <ipc/broad_phase/hash_grid.hpp>
#include <ipc/broad_phase/spatial_hash.hpp>
#include <ipc/broad_phase/brute_force.hpp>

using namespace ipc;

static const bool BENCHMARK_EV = false;
static const bool BENCHMARK_EE = false;
static const bool BENCHMARK_FV = true;

TEST_CASE(
    "Benchmark different spatial hashes",
    "[!benchmark][spatial_hash][hash_grid]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

#ifdef NDEBUG
    std::string filename =
        GENERATE(std::string("cube.obj"), std::string("bunny.obj"));
#else
    std::string filename = "cube.obj";
#endif
    std::string mesh_path = std::string(TEST_DATA_DIR) + filename;
    bool success = igl::read_triangle_mesh(mesh_path, V0, F);
    REQUIRE(success);

    VectorMax3d min_V = V0.colwise().minCoeff();
    VectorMax3d max_V = V0.colwise().maxCoeff();
    VectorMax3d center = (max_V + min_V) / 2;
    V0.rowwise() -= center.transpose();

    SECTION("Squish")
    {
        igl::edges(F, E);
        double s = GENERATE(0.1, 0.0, -0.1);
        V1 = V0;
        V1.col(0) *= s;
    }
    SECTION("Interobject")
    {
        VectorMax3d scale = max_V - min_V;
        for (int j = 0; j < V0.cols(); j++) {
            V0.col(j).array() /= scale(j);
        }

        int n = V0.rows();

        Eigen::MatrixXd V(2 * n, V0.cols());
        V.topRows(n) = V0;
        V.block(0, 0, n, 1).array() += 0.5;
        V.bottomRows(n) = V0;
        V.block(n, 0, n, 1).array() -= 0.5;
        V0 = V;

        double dx = GENERATE(0.0, 0.25, 0.75, 1.0);
        V1 = V0;
        V1.block(0, 0, n, 1).array() -= dx;
        V1.block(n, 0, n, 1).array() += dx;

        Eigen::MatrixXi F_stack(2 * F.rows(), F.cols());
        F_stack.topRows(F.rows()) = F;
        F_stack.bottomRows(F.rows()) = F.array() + n;
        F = F_stack;

        igl::edges(F, E);
    }

    double inflation_radius = 1e-2; // GENERATE(take(5, random(0.0, 0.1)));

    BENCHMARK("SpatialHash")
    {
        SpatialHash sh;
        sh.build(V0, V1, E, F, inflation_radius);

        Candidates candidates;
        sh.queryMeshForCandidates(
            V0, V1, E, F, candidates,
            /*queryEV=*/BENCHMARK_EV,
            /*queryEE=*/BENCHMARK_EE,
            /*queryFV=*/BENCHMARK_FV);
    };

    BENCHMARK("HashGrid")
    {
        HashGrid hashgrid;
        hashgrid.resize(V0, V1, E, inflation_radius);
        hashgrid.addVertices(V0, V1, inflation_radius);
        hashgrid.addEdges(V0, V1, E, inflation_radius);
        hashgrid.addFaces(V0, V1, F, inflation_radius);

        Candidates candidates;
        if (BENCHMARK_EV) {
            hashgrid.getVertexEdgePairs(E, candidates.ev_candidates);
        }
        if (BENCHMARK_EE) {
            hashgrid.getEdgeEdgePairs(E, candidates.ee_candidates);
        }
        if (BENCHMARK_FV) {
            hashgrid.getFaceVertexPairs(F, candidates.fv_candidates);
        }
    };

    BENCHMARK("BruteForce")
    {
        Candidates candidates;
        detect_collision_candidates_brute_force(
            V0, V1, E, F, candidates,
            /*detect_edge_vertex=*/BENCHMARK_EV,
            /*detect_edge_edge=*/BENCHMARK_EE,
            /*detect_face_vertex=*/BENCHMARK_FV,
            /*perform_aabb_check=*/true, inflation_radius);
    };
}

TEST_CASE(
    "Benchmark different spatial hashes on real data",
    "[!benchmark][spatial_hash][hash_grid][real_data]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    std::string mesh_path =
        std::string(TEST_DATA_DIR) + "slow-broadphase-ccd/s0.obj";
    bool success = igl::read_triangle_mesh(mesh_path, V0, F);
    if (!success) {
        return; // Data is private
    }

    mesh_path = std::string(TEST_DATA_DIR) + "slow-broadphase-ccd/s1.obj";
    success = igl::read_triangle_mesh(mesh_path, V1, F);
    if (!success) {
        return; // Data is private
    }

    igl::edges(F, E);

    double inflation_radius = 1e-2; // GENERATE(take(5, random(0.0, 0.1)));

    BENCHMARK("SpatialHash")
    {
        SpatialHash sh;
        sh.build(V0, V1, E, F, inflation_radius);

        Candidates candidates;
        sh.queryMeshForCandidates(
            V0, V1, E, F, candidates,
            /*queryEV=*/BENCHMARK_EV,
            /*queryEE=*/BENCHMARK_EE,
            /*queryFV=*/BENCHMARK_FV);
    };

    BENCHMARK("HashGrid")
    {
        HashGrid hashgrid;
        hashgrid.resize(V0, V1, E, inflation_radius);
        hashgrid.addVertices(V0, V1, inflation_radius);
        hashgrid.addEdges(V0, V1, E, inflation_radius);
        hashgrid.addFaces(V0, V1, F, inflation_radius);

        Candidates candidates;
        if (BENCHMARK_EV) {
            hashgrid.getVertexEdgePairs(E, candidates.ev_candidates);
        }
        if (BENCHMARK_EE) {
            hashgrid.getEdgeEdgePairs(E, candidates.ee_candidates);
        }
        if (BENCHMARK_FV) {
            hashgrid.getFaceVertexPairs(F, candidates.fv_candidates);
        }
    };

    BENCHMARK("BruteForce")
    {
        Candidates candidates;
        detect_collision_candidates_brute_force(
            V0, V1, E, F, candidates,
            /*detect_edge_vertex=*/BENCHMARK_EV,
            /*detect_edge_edge=*/BENCHMARK_EE,
            /*detect_face_vertex=*/BENCHMARK_FV,
            /*perform_aabb_check=*/true, inflation_radius);
    };
}
