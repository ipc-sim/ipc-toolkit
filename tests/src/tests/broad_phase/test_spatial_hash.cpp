#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>

#include <ipc/broad_phase/spatial_hash.hpp>

using namespace ipc;

TEST_CASE("Build SpatialHash", "[broad_phase][spatial_hash][build]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    REQUIRE(tests::load_mesh("cloth_ball92.ply", V0, E, F));
    REQUIRE(tests::load_mesh("cloth_ball93.ply", V1, E, F));

    double inflation_radius = 0;

    SpatialHash sh;
    sh.build(V0, V1, E, F, inflation_radius);

    Candidates candidates;
    sh.detect_collision_candidates(V0.cols(), candidates);

    CHECK(candidates.size() == 6'852'873);
    CHECK(candidates.vv_candidates.size() == 0);
    CHECK(candidates.ev_candidates.size() == 0);
    CHECK(candidates.ee_candidates.size() == 5'197'332);
    CHECK(candidates.fv_candidates.size() == 1'655'541);

    sh.clear();
}
