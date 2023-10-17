#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/config.hpp>
#include <ipc/ccd/ccd.hpp>
#include <ipc/ccd/additive_ccd.hpp>

using namespace ipc;

TEST_CASE("Point-point CCD", "[ccd][point-point]")
{
    const int dim = GENERATE(2, 3);
    const Eigen::Vector3d p0_t0_3D(0.0, 0.0, 0.0), p0_t1_3D(1.0, 1.0, 1.0);
    const Eigen::Vector3d p1_t0_3D(1.0, 1.0, 0.0), p1_t1_3D(0.0, 0.0, 1.0);
    const VectorMax3d p0_t0 = p0_t0_3D.head(dim), p0_t1 = p0_t1_3D.head(dim);
    const VectorMax3d p1_t0 = p1_t0_3D.head(dim), p1_t1 = p1_t1_3D.head(dim);

    const double min_distance = GENERATE(0, 1e-6, 1e-4, 1e-2);

    double toi;
    bool is_colliding =
        point_point_ccd(p0_t0, p1_t0, p0_t1, p1_t1, toi, min_distance);

    // Check the results
    CHECK(is_colliding);
    CHECK(toi == Catch::Approx(0.5).margin(1e-3));

    is_colliding = additive_ccd::point_point_ccd(
        p0_t0, p1_t0, p0_t1, p1_t1, toi, min_distance, /*tmax=*/1,
        /*conservative_rescaling=*/0.999);

    // Check the results
    CHECK(is_colliding);
    CHECK(toi == Catch::Approx(0.5).margin(1e-3));
}