#include <catch2/catch_all.hpp>

#include <Eigen/Core>

#include <igl/IO>
#include <igl/edges.h>

#include <ipc/broad_phase/spatial_hash.hpp>
#include <ipc/utils/logger.hpp>

#include <utils.hpp>

using namespace ipc;

TEST_CASE("Build SpatialHash", "[spatial_hash][build]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    bool success =
        igl::read_triangle_mesh(TEST_DATA_DIR + "cloth_ball92.ply", V0, F);
    REQUIRE(success);

    success =
        igl::read_triangle_mesh(TEST_DATA_DIR + "cloth_ball93.ply", V1, F);
    REQUIRE(success);

    igl::edges(F, E);

    double inflation_radius = 0;

    SpatialHash sh;
    sh.build(V0, V1, E, F, inflation_radius);

    Candidates candidates;
    sh.detect_collision_candidates(V0.cols(), candidates);

    sh.clear();
}
