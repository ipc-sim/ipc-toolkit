#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/candidates/candidates.hpp>

using namespace ipc;

TEST_CASE("Per-vertex safe distances", "[ogc][candidates][safe-distances]")
{
    // Build a simple 3-vertex mesh where vertex 0 and 1 are close and vertex 2
    // is isolated.
    Eigen::MatrixXd V(3, 3);
    V << -1, 0, 0, /**/ 1, 0, 0, /**/ 10, 0, 0;

    Eigen::MatrixXi E, F;

    // One vertex-vertex candidate between 0 and 1
    Candidates C;
    C.vv_candidates = {
        VertexVertexCandidate(0, 1),
        VertexVertexCandidate(0, 2),
        VertexVertexCandidate(1, 2),
    };

    const double inflation_radius = 10.0;
    const double min_distance = 0.5;

    Eigen::VectorXd per_vertex = C.compute_per_vertex_safe_distances(
        /*mesh=*/CollisionMesh(V), V, inflation_radius, min_distance);

    REQUIRE(per_vertex.size() == 3);

    // Distance between vertex 0 and 1 is 2.0 -> safe distance = 2.0 -
    // min_distance
    Eigen::VectorXd expected_per_vertex(per_vertex.size());
    expected_per_vertex(0) = 2 - min_distance;
    expected_per_vertex(1) = 2 - min_distance;
    expected_per_vertex(2) = 9 - min_distance;

    CHECK(per_vertex[0] == Catch::Approx(expected_per_vertex[0]));
    CHECK(per_vertex[1] == Catch::Approx(expected_per_vertex[1]));
    CHECK(per_vertex[2] == Catch::Approx(expected_per_vertex[2]));

    // All entries must be non-negative
    for (int i = 0; i < per_vertex.size(); ++i) {
        CHECK(per_vertex[i] >= 0.0);
    }
}
