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
    "Nearest rotation of a rectangular matrix", "[math][nearest_rotation]")
{
    // A tall matrix has orthonormal columns and a wide one has orthonormal
    // rows. Either way, we build it as an orthonormal frame times a symmetric
    // positive definite stretch, so the frame is the nearest rotation.
    const bool tall = GENERATE(true, false);
    const int rows = tall ? 3 : 2;
    const int cols = tall ? 2 : 3;

    const Eigen::Matrix3d R =
        Eigen::AngleAxisd(0.7, Eigen::Vector3d(1, 2, 3).normalized())
            .toRotationMatrix();
    const MatrixMax3d Q =
        tall ? MatrixMax3d(R.leftCols(2)) : MatrixMax3d(R.topRows(2));

    MatrixMax3d B = MatrixMax3d::Zero(2, 2);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            B(i, j) = std::sin(1.0 + 3 * i + j);
        }
    }
    const MatrixMax3d P = MatrixMax3d::Identity(2, 2) + 0.5 * B.transpose() * B;
    const MatrixMax3d A = tall ? MatrixMax3d(Q * P) : MatrixMax3d(P * Q);

    SECTION("Recovers the frame of a stretched frame")
    {
        const MatrixMax3d N = nearest_rotation(A);
        CHECK(N.rows() == rows);
        CHECK(N.cols() == cols);
        CHECK((N - Q).norm() < 1e-12);
    }

    SECTION("Is unaffected by a reflection of the input's smaller side")
    {
        // Flipping one axis of the 2D side changes the orientation of the
        // frame, but the result must still be a valid orthonormal frame.
        MatrixMax3d reflection = MatrixMax3d::Identity(2, 2);
        reflection(1, 1) = -1;
        const MatrixMax3d flipped =
            tall ? MatrixMax3d(A * reflection) : MatrixMax3d(reflection * A);

        const MatrixMax3d N = nearest_rotation(flipped);
        const MatrixMax3d gram = tall ? MatrixMax3d(N.transpose() * N)
                                      : MatrixMax3d(N * N.transpose());
        CHECK((gram - MatrixMax3d::Identity(2, 2)).norm() < 1e-12);
        CHECK(
            (N
             - (tall ? MatrixMax3d(Q * reflection)
                     : MatrixMax3d(reflection * Q)))
                .norm()
            < 1e-12);
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
