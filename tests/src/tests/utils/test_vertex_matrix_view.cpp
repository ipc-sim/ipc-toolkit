#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/esp/collisions/vertex_matrix_view.hpp>

using namespace ipc;
using Catch::Approx;

TEST_CASE("VertexMatrixView single matrix", "[vertex_matrix_view]")
{
    Eigen::MatrixXd A(3, 3);
    // clang-format off
    A << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
    // clang-format on

    VertexMatrixView<3> view(A);

    REQUIRE(view.rows() == 3);
    REQUIRE(view.cols() == 3);
    REQUIRE(view.m_B == nullptr);

    for (index_t i = 0; i < 3; i++) {
        auto row = view(i);
        for (int j = 0; j < 3; j++) {
            CHECK(row(j) == Approx(A(i, j)));
        }
    }
}

TEST_CASE("VertexMatrixView two-matrix concatenation", "[vertex_matrix_view]")
{
    Eigen::MatrixXd A(2, 3);
    Eigen::MatrixXd B(3, 3);
    // clang-format off
    A << 1, 2, 3,
         4, 5, 6;
    B << 7,  8,  9,
         10, 11, 12,
         13, 14, 15;
    // clang-format on

    VertexMatrixView<3> view(A, B);

    REQUIRE(view.rows() == 5);
    REQUIRE(view.cols() == 3);
    REQUIRE(view.n_A_rows == 2);
    REQUIRE(view.n_B_rows == 3);

    // Check rows from A
    for (index_t i = 0; i < 2; i++) {
        auto row = view(i);
        for (int j = 0; j < 3; j++) {
            CHECK(row(j) == Approx(A(i, j)));
        }
    }

    // Check rows from B
    for (index_t i = 0; i < 3; i++) {
        auto row = view(2 + i);
        for (int j = 0; j < 3; j++) {
            CHECK(row(j) == Approx(B(i, j)));
        }
    }
}

TEST_CASE(
    "VertexMatrixView concatenation matches naive vstack",
    "[vertex_matrix_view]")
{
    // Build random matrices and verify the view matches manual stacking
    const int nA = 4, nB = 5;
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(nA, 3);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(nB, 3);

    Eigen::MatrixXd AB(nA + nB, 3);
    AB.topRows(nA) = A;
    AB.bottomRows(nB) = B;

    VertexMatrixView<3> view(A, B);

    REQUIRE(view.rows() == nA + nB);

    for (index_t i = 0; i < view.rows(); i++) {
        auto row = view(i);
        for (int j = 0; j < 3; j++) {
            CHECK(row(j) == Approx(AB(i, j)));
        }
    }
}

TEST_CASE("VertexMatrixView non-owning semantics", "[vertex_matrix_view]")
{
    Eigen::MatrixXd A(2, 3);
    Eigen::MatrixXd B(1, 3);
    A << 1, 2, 3, 4, 5, 6;
    B << 7, 8, 9;

    VertexMatrixView<3> view(A, B);

    // The view should point to the original data
    CHECK(view.m_A == A.data());
    CHECK(view.m_B == B.data());

    // Mutate A and verify the view reflects the change
    A(0, 0) = 99;
    auto row = view(0);
    CHECK(row(0) == Approx(99));
}

TEST_CASE("VertexMatrixView empty matrix B", "[vertex_matrix_view]")
{
    Eigen::MatrixXd A(3, 3);
    A << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    Eigen::MatrixXd B(0, 3);

    VertexMatrixView<3> view(A, B);

    REQUIRE(view.rows() == 3);
    REQUIRE(view.n_B_rows == 0);

    for (index_t i = 0; i < 3; i++) {
        auto row = view(i);
        for (int j = 0; j < 3; j++) {
            CHECK(row(j) == Approx(A(i, j)));
        }
    }
}
