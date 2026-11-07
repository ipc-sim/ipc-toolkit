#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/geometry/normal.hpp>
#include <ipc/math/math.hpp>

#include <Eigen/Geometry>

using namespace ipc;

TEST_CASE("Nearest rotation of a matrix", "[math][nearest_rotation]")
{
    const int dim = GENERATE(2, 3);

    MatrixMax3d R;
    if (dim == 2) {
        R = Eigen::Rotation2Dd(0.7).toRotationMatrix();
    } else {
        R = Eigen::AngleAxisd(0.7, Eigen::Vector3d(1, 2, 3).normalized())
                .toRotationMatrix();
    }

    // The nearest rotation to the polar decomposition A = R P (with P
    // symmetric positive definite) is exactly its rotation factor R.
    MatrixMax3d B = MatrixMax3d::Zero(dim, dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            B(i, j) = std::sin(1.0 + 3 * i + j);
        }
    }
    const MatrixMax3d P =
        MatrixMax3d::Identity(dim, dim) + 0.5 * B.transpose() * B;

    SECTION("Recovers the rotation of a stretched rotation")
    {
        CHECK((nearest_rotation(R * P) - R).norm() < 1e-12);
    }

    SECTION("Always returns a proper rotation, even for a reflection")
    {
        MatrixMax3d reflection = MatrixMax3d::Identity(dim, dim);
        reflection(dim - 1, dim - 1) = -1;

        const MatrixMax3d Q = nearest_rotation(R * P * reflection);
        CHECK(
            (Q.transpose() * Q - MatrixMax3d::Identity(dim, dim)).norm()
            < 1e-12);
        CHECK(Q.determinant() == Catch::Approx(1));
    }
}

TEST_CASE(
    "2D cross product matrix acts as a scalar angular velocity",
    "[geometry][cross_product_matrix]")
{
    // In 2D, ω × v is ω times v rotated a quarter turn counter-clockwise.
    const double omega = 0.5;
    const Eigen::Vector2d v(2, 3);
    CHECK(
        (cross_product_matrix(omega) * v
         - omega * Eigen::Vector2d(-v.y(), v.x()))
            .norm()
        < 1e-15);
}
