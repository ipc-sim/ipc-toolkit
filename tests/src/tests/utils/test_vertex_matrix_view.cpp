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
    REQUIRE(view.m_b == nullptr);

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
    REQUIRE(view.m_n_a_rows == 2);
    REQUIRE(view.m_n_b_rows == 3);

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
    CHECK(view.m_a == A.data());
    CHECK(view.m_b == B.data());

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
    REQUIRE(view.m_n_b_rows == 0);

    for (index_t i = 0; i < 3; i++) {
        auto row = view(i);
        for (int j = 0; j < 3; j++) {
            CHECK(row(j) == Approx(A(i, j)));
        }
    }
}

TEST_CASE(
    "VertexMatrixView reads strided and row-major matrices",
    "[vertex_matrix_view]")
{
    const Eigen::MatrixXd A = Eigen::MatrixXd::Random(4, 3);
    const Eigen::MatrixXd B = Eigen::MatrixXd::Random(2, 3);
    Eigen::MatrixXd AB(6, 3);
    AB << A, B;

    // Embed A and B in taller matrices so each column is further from the
    // next than the number of rows the view sees.
    Eigen::MatrixXd A_big = Eigen::MatrixXd::Random(10, 3);
    A_big.topRows(4) = A;
    Eigen::MatrixXd B_big = Eigen::MatrixXd::Random(3, 5);
    B_big.block<2, 3>(1, 2) = B;
    const Eigen::Matrix<double, -1, 3, Eigen::RowMajor> A_row_major = A;

    const auto check_matches = [](const VertexMatrixView<3>& view,
                                  const Eigen::MatrixXd& expected) {
        REQUIRE(view.rows() == expected.rows());
        for (index_t i = 0; i < view.rows(); i++) {
            CHECK(view(i) == expected.row(i));
        }
    };

    SECTION("Block")
    {
        check_matches(VertexMatrixView<3>(A_big.topRows(4)), A);
        check_matches(
            VertexMatrixView<3>(A_big.topRows(4), B_big.block<2, 3>(1, 2)), AB);
    }
    SECTION("Ref")
    {
        // A const Ref binds blocks without a copy, so the view must use the
        // Ref's outer stride.
        const Eigen::Ref<const Eigen::MatrixXd> A_ref = A_big.topRows(4);
        REQUIRE(A_ref.data() == A_big.data());
        const VertexMatrixView<3> view(A_ref, B);
        CHECK(view.m_a == A_big.data());
        check_matches(view, AB);
    }
    SECTION("Map with outer stride")
    {
        const Eigen::Map<const Eigen::MatrixXd, 0, Eigen::OuterStride<>> A_map(
            A_big.data(), 4, 3, Eigen::OuterStride<>(A_big.rows()));
        check_matches(VertexMatrixView<3>(A_map), A);
    }
    SECTION("Row-major")
    {
        const VertexMatrixView<3> view(A_row_major, B);
        CHECK(view.m_a == A_row_major.data());
        check_matches(view, AB);
    }
    SECTION("Row vector")
    {
        const Eigen::RowVector3d q = B.row(1);
        const VertexMatrixView<3> view(A_big.topRows(4), q);
        CHECK(view.m_b == q.data());
        Eigen::MatrixXd Aq(5, 3);
        Aq << A, q;
        check_matches(view, Aq);
    }
}
