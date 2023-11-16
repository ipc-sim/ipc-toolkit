#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/broad_phase/aabb.hpp>

using namespace ipc;

// TEST_CASE("AABB initilization", "[broad_phase][AABB]")
// {
//     int dim = GENERATE(2, 3);
//     CAPTURE(dim);
//     AABB aabb;
//     ArrayMax3d actual_center(dim);
//     SECTION("Empty AABB")
//     {
//         aabb = AABB(ArrayMax3d::Zero(dim), ArrayMax3d::Zero(dim));
//         actual_center.setZero();
//     }
//     SECTION("Box centered at zero")
//     {
//         ArrayMax3d min =
//             ArrayMax3d::Random(dim).array() - 1; // in range [-2, 0]
//         ArrayMax3d max = -min;
//         aabb = AABB(min, max);
//         actual_center.setZero();
//     }
//     SECTION("Box not centered at zero")
//     {
//         ArrayMax3d min(dim), max(dim);
//         if (dim == 2) {
//             min << 5.1, 3.14;
//             max << 10.4, 7.89;
//             actual_center << 7.75, 5.515;
//         } else {
//             min << 5.1, 3.14, 7.94;
//             max << 10.4, 7.89, 10.89;
//             actual_center << 7.75, 5.515, 9.415;
//         }
//         aabb = AABB(min, max);
//     }
//     ArrayMax3d center_diff = aabb.getCenter() - actual_center;
//     CHECK(center_diff.matrix().norm() == Catch::Approx(0.0).margin(1e-12));
// }

TEST_CASE("AABB overlapping", "[broad_phase][AABB]")
{
    AABB a, b;
    bool are_overlapping = false;
    SECTION("a to the right of b")
    {
        a = AABB(Eigen::Array2d(-1, 0), Eigen::Array2d(0, 1));
        SECTION("overlapping")
        {
            b = AABB(Eigen::Array2d(-0.5, 0), Eigen::Array2d(0.5, 1));
            are_overlapping = true;
        }
        SECTION("not overlapping")
        {
            b = AABB(Eigen::Array2d(0.5, 0), Eigen::Array2d(1.5, 1));
            are_overlapping = false;
        }
    }
    SECTION("b to the right of a")
    {
        b = AABB(Eigen::Array2d(-1, 0), Eigen::Array2d(0, 1));
        SECTION("overlapping")
        {
            a = AABB(Eigen::Array2d(-0.5, 0), Eigen::Array2d(0.5, 1));
            are_overlapping = true;
        }
        SECTION("not overlapping")
        {
            a = AABB(Eigen::Array2d(0.5, 0), Eigen::Array2d(1.5, 1));
            are_overlapping = false;
        }
    }
    SECTION("a above b")
    {
        a = AABB(Eigen::Array2d(0, -1), Eigen::Array2d(1, 0));
        SECTION("overlapping")
        {
            b = AABB(Eigen::Array2d(0, -0.5), Eigen::Array2d(1, 0.5));
            are_overlapping = true;
        }
        SECTION("not overlapping")
        {
            b = AABB(Eigen::Array2d(0, 0.5), Eigen::Array2d(1, 1.5));
            are_overlapping = false;
        }
    }
    SECTION("a above b")
    {
        b = AABB(Eigen::Array2d(0, -1), Eigen::Array2d(1, 0));
        SECTION("overlapping")
        {
            a = AABB(Eigen::Array2d(0, -0.5), Eigen::Array2d(1, 0.5));
            are_overlapping = true;
        }
        SECTION("not overlapping")
        {
            a = AABB(Eigen::Array2d(0, 0.5), Eigen::Array2d(1, 1.5));
            are_overlapping = false;
        }
    }
    CHECK(a.intersects(b) == are_overlapping);
}
