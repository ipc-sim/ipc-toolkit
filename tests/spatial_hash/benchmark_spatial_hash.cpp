#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <igl/IO>
#include <igl/edges.h>
#include <igl/Timer.h>

#include <ipc/ipc.hpp>
#include <ipc/spatial_hash/hash_grid.hpp>
#include <ipc/spatial_hash/spatial_hash.hpp>
#include <ipc/spatial_hash/brute_force.hpp>

using namespace ipc;

TEST_CASE(
    "Benchmark different spatial hashes",
    "[!benchmark][spatial_hash][hash_grid]")
{
    using namespace ipc;

    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;
    Eigen::VectorXi group_ids;

#ifdef NDEBUG
    std::string filename =
        GENERATE(std::string("cube.obj"), std::string("bunny.obj"));
#else
    std::string filename = "cube.obj";
#endif
    std::string mesh_path = std::string(TEST_DATA_DIR) + filename;
    bool success = igl::read_triangle_mesh(mesh_path, V0, F);
    REQUIRE(success);

    Eigen::VectorX3d min_V = V0.colwise().minCoeff();
    Eigen::VectorX3d max_V = V0.colwise().maxCoeff();
    Eigen::VectorX3d center = (max_V + min_V) / 2;
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
        Eigen::VectorX3d scale = max_V - min_V;
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
            /*queryEV=*/true, /*queryEE=*/true, /*queryFV=*/true);
    };

    BENCHMARK("HashGrid")
    {
        HashGrid hashgrid;
        hashgrid.resize(V0, V1, E, inflation_radius);
        hashgrid.addVertices(V0, V1, inflation_radius);
        hashgrid.addEdges(V0, V1, E, inflation_radius);
        hashgrid.addFaces(V0, V1, F, inflation_radius);

        Candidates candidates;
        hashgrid.getVertexEdgePairs(E, group_ids, candidates.ev_candidates);
        hashgrid.getEdgeEdgePairs(E, group_ids, candidates.ee_candidates);
        hashgrid.getFaceVertexPairs(F, group_ids, candidates.fv_candidates);
    };

    BENCHMARK("BruteForce")
    {
        Candidates candidates;
        detect_collision_candidates_brute_force(
            V0, V1, E, F, candidates, /*detect_edge_vertex=*/true,
            /*detect_edge_edge=*/true, /*detect_face_vertex=*/true,
            /*perform_aabb_check=*/true, inflation_radius, group_ids);
    };
}
