#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <igl/IO>
#include <igl/edges.h>
#include <igl/Timer.h>

#include <ipc/ipc.hpp>
#include <ipc/spatial_hash/hash_grid.hpp>
#include <ipc/spatial_hash/spatial_hash.hpp>

using namespace ipc;

TEST_CASE(
    "Benchmark different spatial hashes",
    "[!benchmark][spatial_hash][hash_grid]")
{
    using namespace ipc;

    Eigen::MatrixXd V, U;
    Eigen::MatrixXi E, F;
    Eigen::VectorXi group_ids;

    SECTION("Simple")
    {
        V.resize(4, 3);
        V.row(0) << -1, -1, 0;
        V.row(1) << 1, -1, 0;
        V.row(2) << 0, 1, 1;
        V.row(3) << 0, 1, -1;

        E.resize(2, 2);
        E.row(0) << 0, 1;
        E.row(1) << 2, 3;

        SECTION("Without group ids") {}
        SECTION("With group ids")
        {
            group_ids.resize(4);
            group_ids << 0, 0, 1, 1;
        }

        F.resize(0, 3);

        U = Eigen::MatrixXd::Zero(V.rows(), V.cols());
        U.col(1).head(2).setConstant(2);
        U.col(1).tail(2).setConstant(-2);
    }
    SECTION("Complex")
    {
#ifdef NDEBUG
        std::string filename =
            GENERATE(std::string("cube.obj"), std::string("bunny.obj"));
#else
        std::string filename = "cube.obj";
#endif
        std::string mesh_path = std::string(TEST_DATA_DIR) + filename;
        bool success = igl::read_triangle_mesh(mesh_path, V, F);
        REQUIRE(success);
        igl::edges(F, E);

        U = Eigen::MatrixXd::Zero(V.rows(), V.cols());
        U.col(1).setOnes();
    }

    HashGrid hashgrid;
    SpatialHash sh;
    Candidates candidates;

    double inflation_radius = 1e-2; // GENERATE(take(5, random(0.0, 0.1)));

    for (int i = 0; i < 2; i++) {
        BENCHMARK("SpatialHash")
        {
            sh.build(V, V + U, E, F);

            candidates.clear();

            sh.queryMeshForCandidates(
                V, V + U, E, F, candidates, inflation_radius,
                /*queryEV=*/true, /*queryEE=*/true, /*queryFV=*/true);
        };
        BENCHMARK("HashGrid")
        {
            hashgrid.resize(V, V + U, E, inflation_radius);
            hashgrid.addVertices(V, V + U, inflation_radius);
            hashgrid.addEdges(V, V + U, E, inflation_radius);
            hashgrid.addFaces(V, V + U, F, inflation_radius);

            candidates.clear();

            hashgrid.getVertexEdgePairs(E, group_ids, candidates.ev_candidates);
            hashgrid.getEdgeEdgePairs(E, group_ids, candidates.ee_candidates);
            hashgrid.getFaceVertexPairs(F, group_ids, candidates.fv_candidates);
        };

        U.setRandom();
        U *= 3;
    }
}
