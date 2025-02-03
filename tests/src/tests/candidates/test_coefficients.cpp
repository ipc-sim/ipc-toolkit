#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/candidates/candidates.hpp>

#include <igl/edges.h>

using namespace ipc;

TEST_CASE("Face-Vertex collision stencil coeffs.", "[fv][stencil][coeffs]")
{
    Eigen::MatrixXd V(4, 3);
    V << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1;

    SECTION("default") { }
    SECTION("random") { V.row(3).setRandom(); }
    SECTION("off triangle 0") { V.row(3) << 2, 0.5, 1; }
    SECTION("off triangle 1") { V.row(3) << 2, 2, 1; }
    SECTION("off triangle 2") { V.row(3) << 0.5, 2, 1; }
    SECTION("off triangle 3") { V.row(3) << -1, 2, 1; }
    SECTION("off triangle 4") { V.row(3) << -1, 0.5, 1; }
    SECTION("off triangle 5") { V.row(3) << -1, -1, 1; }
    SECTION("off triangle 6") { V.row(3) << 0.5, -1, 1; }
    SECTION("off triangle 7") { V.row(3) << 2, -0.5, 1; }

    Eigen::MatrixXi F(1, 3);
    F << 0, 1, 2;

    Eigen::MatrixXi E;
    igl::edges(F, E);

    FaceVertexCandidate fv(0, 3);

    VectorMax4d coeffs = fv.compute_coefficients(V, E, F);
    CAPTURE(V.row(3), coeffs.transpose());

    CHECK(-coeffs[1] == Catch::Approx(1 + coeffs[2] + coeffs[3]));

    std::array<VectorMax3d, 4> vertices = fv.vertices(V, E, F);
    Eigen::Vector3d n = Eigen::Vector3d::Zero();
    for (int i = 0; i < fv.num_vertices(); i++) {
        n += coeffs[i] * vertices[i];
    }

    CHECK(n.squaredNorm() == Catch::Approx(fv.compute_distance(V, E, F)));
}

TEST_CASE("Edge-edge collision stencil coeffs.", "[ee][stencil][coeffs]")
{
    Eigen::MatrixXd V(4, 3);
    V << -1, 0, 0, /**/ 1, 0, 0, /**/ 0, -1, 0, /**/ 0, 1, 0;

    // clang-format off
    SECTION("default") { }
    SECTION("random") { V.bottomRows<2>().setRandom(); }
    SECTION("quad 0") { V.block<2, 1>(2, 0).array() += 2; }
    SECTION("quad 1") { V.block<2, 2>(2, 0).array() += 2; }
    SECTION("quad 2") { V.block<2, 1>(2, 1).array() += 2; }
    SECTION("quad 3") { V.bottomRows<2>().rowwise() += Eigen::RowVector3d(-2, 2, 0); }
    SECTION("quad 4") { V.block<2, 1>(2, 0).array() -= 2; }
    SECTION("quad 5") { V.block<2, 2>(2, 0).array() -= 2; }
    SECTION("quad 6") { V.block<2, 1>(2, 1).array() -= 2; }
    SECTION("quad 7") { V.bottomRows<2>().rowwise() += Eigen::RowVector3d(2, -2, 0); }
    // clang-format on

    Eigen::MatrixXi E(2, 2), F;
    E << 0, 1, /**/ 2, 3;

    EdgeEdgeCandidate ee(0, 1);

    VectorMax4d coeffs = ee.compute_coefficients(V, E, F);
    CAPTURE(
        V, coeffs.transpose(),
        edge_edge_distance_type(V.row(0), V.row(1), V.row(2), V.row(3)));

    std::array<VectorMax3d, 4> vertices = ee.vertices(V, E, F);
    Eigen::Vector3d n = Eigen::Vector3d::Zero();
    for (int i = 0; i < ee.num_vertices(); i++) {
        n += coeffs[i] * vertices[i];
    }

    CHECK(n.squaredNorm() == Catch::Approx(ee.compute_distance(V, E, F)));
}