#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/candidates/candidates.hpp>
#include <ipc/geometry/normal.hpp>
#include <ipc/collisions/normal/plane_vertex.hpp>

#include <igl/edges.h>

#include <finitediff.hpp>

using namespace ipc;

TEST_CASE("Vertex-vertex collision normal", "[vv][normal]")
{
    const int dim = GENERATE(2, 3);

    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(2, dim);
    V(0, 0) = -1;
    V(1, 0) = 1;

    VectorMax3d expected_normal = VectorMax3d::Zero(dim);
    expected_normal(0) = -1;

    Eigen::MatrixXi E, F;

    VertexVertexCandidate vv(0, 1);

    VectorMax3d normal = vv.compute_normal(V, E, F);
    MatrixMax<double, 3, 12> jacobian = vv.compute_normal_jacobian(V, E, F);

    REQUIRE(normal.size() == dim);
    REQUIRE(jacobian.rows() == dim);
    REQUIRE(jacobian.cols() == 2 * dim);

    CHECK(normal.norm() == Catch::Approx(1.0));
    CHECK(normal.isApprox(expected_normal));

    // Check jacobian using finite differences
    Eigen::MatrixXd fd_jacobian;
    VectorMax12d x = vv.dof(V, E, F);
    fd::finite_jacobian(
        x,
        [&vv](const Eigen::MatrixXd& x_fd) { return vv.compute_normal(x_fd); },
        fd_jacobian);
    CHECK(fd::compare_jacobian(jacobian, fd_jacobian));
}

TEST_CASE("Edge-vertex collision normal", "[ev][normal]")
{
    const int dim = GENERATE(2, 3);

    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(3, dim);
    V(0, 0) = -1; // Vertex 0 x=-1
    V(1, 0) = 1;  // Vertex 1 x=1
    V(2, 1) = 1;  // Vertex 2 y=1

    SECTION("default") { }
    SECTION("flip") { V(2, 1) = -1; }
    SECTION("e0") { V(2, 0) = -2; }
    SECTION("e1") { V(2, 0) = 2; }

    CAPTURE(dim, V(2, 0), V(2, 1));

    Eigen::MatrixXi E(1, 2), F;
    E << 0, 1;

    EdgeVertexCandidate ev(0, 2);

    auto normal = ev.compute_normal(V, E, F);
    auto jacobian = ev.compute_normal_jacobian(V, E, F);

    REQUIRE(normal.size() == dim);
    REQUIRE(jacobian.rows() == dim);
    REQUIRE(jacobian.cols() == 3 * dim);

    CHECK(normal.norm() == Catch::Approx(1.0));
    VectorMax3d expected_normal = VectorMax3d::Zero(dim);
    expected_normal(1) = V(2, 1) < 0 ? -1 : 1;
    CHECK(normal.isApprox(expected_normal));
    if (!normal.isApprox(expected_normal)) {
        std::cout << "Normal: " << normal.transpose() << std::endl;
        std::cout << "Expected: " << expected_normal.transpose() << std::endl;
    }

    // Check jacobian using finite differences
    Eigen::MatrixXd fd_jacobian;
    VectorMax12d x = ev.dof(V, E, F);
    fd::finite_jacobian(
        x,
        [&ev](const Eigen::MatrixXd& x_fd) { return ev.compute_normal(x_fd); },
        fd_jacobian);
    CHECK(fd::compare_jacobian(jacobian, fd_jacobian));
    if (!fd::compare_jacobian(jacobian, fd_jacobian)) {
        std::cout << "Jacobian:\n" << jacobian << std::endl;
        std::cout << "FD Jacobian:\n" << fd_jacobian << std::endl;
    }
}

TEST_CASE("Point-line normal hessian", "[pl][normal]")
{
    const int DIM = GENERATE(2, 3);

    VectorMax3d p(DIM);
    VectorMax3d e0(DIM);
    VectorMax3d e1(DIM);
    Eigen::VectorXd x(3 * DIM);

    const int case_i = GENERATE(range(0, 10));
    p.setRandom();
    e0.setRandom();
    e1.setRandom();

    x << p, e0, e1;

    // Check hessian using finite differences
    Eigen::MatrixXd hessian = point_line_unnormalized_normal_hessian(p, e0, e1);
    Eigen::MatrixXd fd_hessian;
    fd::finite_jacobian(
        x,
        [DIM](const Eigen::VectorXd& x_fd) -> Eigen::MatrixXd {
            return point_line_unnormalized_normal_jacobian(
                x_fd.segment(0, DIM), x_fd.segment(DIM, DIM),
                x_fd.segment(2 * DIM, DIM));
        },
        fd_hessian);
    CHECK(fd::compare_jacobian(hessian, fd_hessian, 1e-6));
    if (!fd::compare_jacobian(hessian, fd_hessian, 1e-6)) {
        std::cout << "Hessian:\n" << hessian << std::endl;
        std::cout << "FD Hessian:\n" << fd_hessian << std::endl;
    }

    // Check hessian using finite differences
    hessian = point_line_normal_hessian(p, e0, e1);
    fd::finite_jacobian(
        x,
        [DIM](const Eigen::VectorXd& x_fd) -> Eigen::MatrixXd {
            return point_line_normal_jacobian(
                x_fd.segment(0, DIM), x_fd.segment(DIM, DIM),
                x_fd.segment(2 * DIM, DIM));
        },
        fd_hessian);
    CHECK(fd::compare_jacobian(hessian, fd_hessian, 1e-6));
    if (!fd::compare_jacobian(hessian, fd_hessian, 1e-6)) {
        std::cout << "Hessian:\n" << hessian << std::endl;
        std::cout << "FD Hessian:\n" << fd_hessian << std::endl;
    }
}

TEST_CASE("Edge-edge collision normal", "[ee][normal]")
{
    Eigen::MatrixXd V(4, 3);
    V << -1, 0, -0.1, /**/ 1, 0, -0.1, /**/ 0, -1, 0.1, /**/ 0, 1, 0.1;

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

    Eigen::Vector3d normal = ee.compute_normal(V, E, F);
    Eigen::MatrixXd jacobian = ee.compute_normal_jacobian(V, E, F);

    REQUIRE(jacobian.rows() == 3);
    REQUIRE(jacobian.cols() == 12);

    CHECK(normal.norm() == Catch::Approx(1.0));

    // Check jacobian using finite differences
    Eigen::MatrixXd fd_jacobian;
    VectorMax12d x = ee.dof(V, E, F);
    fd::finite_jacobian(
        x,
        [&ee](const Eigen::MatrixXd& x_fd) { return ee.compute_normal(x_fd); },
        fd_jacobian);
    CHECK(fd::compare_jacobian(jacobian, fd_jacobian));
    if (!fd::compare_jacobian(jacobian, fd_jacobian)) {
        std::cout << "Jacobian:\n" << jacobian << std::endl;
        std::cout << "FD Jacobian:\n" << fd_jacobian << std::endl;
    }
}

TEST_CASE("Line-line normal hessian", "[ee][normal][hessian]")
{
    Eigen::Vector3d a(0, 0, 0);
    Eigen::Vector3d b(1, 0, 0);
    Eigen::Vector3d c(0, 1, 0);
    Eigen::Vector3d d(0, 1, 1);
    Eigen::VectorXd x(12);

    const int case_i = GENERATE(range(0, 10));
    if (case_i > 0) {
        a.setRandom();
        b.setRandom();
        c.setRandom();
        d.setRandom();
    }

    x << a, b, c, d;

    // Check hessian using finite differences
    Eigen::MatrixXd hessian = line_line_unnormalized_normal_hessian(a, b, c, d);
    Eigen::MatrixXd fd_hessian;
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& x_fd) -> Eigen::MatrixXd {
            return line_line_unnormalized_normal_jacobian(
                x_fd.segment<3>(0), x_fd.segment<3>(3), x_fd.segment<3>(6),
                x_fd.segment<3>(9));
        },
        fd_hessian);
    CHECK(fd::compare_jacobian(hessian, fd_hessian, 1e-6));
    if (!fd::compare_jacobian(hessian, fd_hessian, 1e-6)) {
        std::cout << "Hessian:\n" << hessian << std::endl;
        std::cout << "FD Hessian:\n" << fd_hessian << std::endl;
    }

    // Check hessian using finite differences
    hessian = line_line_normal_hessian(a, b, c, d);
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& x_fd) -> Eigen::MatrixXd {
            return line_line_normal_jacobian(
                x_fd.segment<3>(0), x_fd.segment<3>(3), x_fd.segment<3>(6),
                x_fd.segment<3>(9));
        },
        fd_hessian);
    CHECK(fd::compare_jacobian(hessian, fd_hessian, 1e-6));
    if (!fd::compare_jacobian(hessian, fd_hessian, 1e-6)) {
        std::cout << "Hessian:\n" << hessian << std::endl;
        std::cout << "FD Hessian:\n" << fd_hessian << std::endl;
    }
}

TEST_CASE("Face-vertex collision normal", "[fv][normal]")
{
    Eigen::MatrixXd V(4, 3);
    V << 0, 0, 0, /**/ 1, 0, 0, /**/ 0, 1, 0, /**/ 0.333, 0.333, 1;

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
    SECTION("failure case 1") { V.row(3) << 0.680375, -0.211234, 0.566198; }
    SECTION("failure case 2") { V.row(3) << -0.997497, 0.127171, -0.613392; }

    Eigen::MatrixXi F(1, 3);
    F << 0, 1, 2;

    Eigen::MatrixXi E;
    igl::edges(F, E);

    FaceVertexCandidate fv(0, 3);

    Eigen::Vector3d normal = fv.compute_normal(V, E, F);
    Eigen::MatrixXd jacobian = fv.compute_normal_jacobian(V, E, F);

    REQUIRE(jacobian.rows() == 3);
    REQUIRE(jacobian.cols() == 12);

    CHECK(normal.norm() == Catch::Approx(1.0));

    // Check jacobian using finite differences
    Eigen::MatrixXd fd_jacobian;
    VectorMax12d x = fv.dof(V, E, F);
    fd::finite_jacobian(
        x,
        [&fv](const Eigen::MatrixXd& x_fd) { return fv.compute_normal(x_fd); },
        fd_jacobian);
    CHECK(fd::compare_jacobian(jacobian, fd_jacobian));
    if (!fd::compare_jacobian(jacobian, fd_jacobian)) {
        std::cout << "Jacobian:\n" << jacobian << std::endl;
        std::cout << "FD Jacobian:\n" << fd_jacobian << std::endl;
    }
}

TEST_CASE("Triangle normal hessian", "[normal]")
{
    Eigen::Vector3d a(0, 0, 0);
    Eigen::Vector3d b(1, 0, 0);
    Eigen::Vector3d c(0, 1, 0);
    Eigen::VectorXd x(9);

    const int case_i = GENERATE(range(0, 10));
    if (case_i > 0) {
        a.setRandom();
        b.setRandom();
        c.setRandom();
    }

    x << a, b, c;

    // Cross product matrix jacobian
    Eigen::MatrixXd J_cross = cross_product_matrix_jacobian();
    Eigen::MatrixXd fd_J_cross;
    fd::finite_jacobian(
        a,
        [](const Eigen::Vector3d& a_fd) { return cross_product_matrix(a_fd); },
        fd_J_cross);
    CHECK(fd::compare_jacobian(J_cross, fd_J_cross, 1e-6));
    if (!fd::compare_jacobian(J_cross, fd_J_cross, 1e-6)) {
        std::cout << "Hessian:\n" << J_cross << std::endl;
        std::cout << "FD Hessian:\n" << fd_J_cross << std::endl;
    }

    // Check hessian using finite differences
    Eigen::MatrixXd hessian = triangle_unnormalized_normal_hessian(a, b, c);
    Eigen::MatrixXd fd_hessian;
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& x_fd) -> Eigen::MatrixXd {
            return triangle_unnormalized_normal_jacobian(
                x_fd.segment<3>(0), x_fd.segment<3>(3), x_fd.segment<3>(6));
        },
        fd_hessian);
    CHECK(fd::compare_jacobian(hessian, fd_hessian, 1e-6));
    if (!fd::compare_jacobian(hessian, fd_hessian, 1e-6)) {
        std::cout << "Hessian:\n" << hessian << std::endl;
        std::cout << "FD Hessian:\n" << fd_hessian << std::endl;
    }

    // Check hessian using finite differences
    hessian = triangle_normal_hessian(a, b, c);
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& x_fd) -> Eigen::MatrixXd {
            return triangle_normal_jacobian(
                x_fd.segment<3>(0), x_fd.segment<3>(3), x_fd.segment<3>(6));
        },
        fd_hessian);
    CHECK(fd::compare_jacobian(hessian, fd_hessian, 1e-6));
    if (!fd::compare_jacobian(hessian, fd_hessian, 1e-6)) {
        std::cout << "Hessian:\n" << hessian << std::endl;
        std::cout << "FD Hessian:\n" << fd_hessian << std::endl;
    }
}

TEST_CASE("Plane-vertex collision normal", "[pv][normal]")
{
    Eigen::MatrixXd V(1, 3);
    V << 0, 1, 0;

    Eigen::Vector3d n(0, 1, 0), o(0, 0, 0);

    Eigen::MatrixXi E, F;

    PlaneVertexNormalCollision pv(o, n, 0);

    Eigen::Vector3d normal = pv.compute_normal(V, E, F);
    Eigen::MatrixXd jacobian = pv.compute_normal_jacobian(V, E, F);

    CHECK(normal.isApprox(n));
    CHECK(jacobian.rows() == 3);
    CHECK(jacobian.cols() == 3);

    // Check jacobian using finite differences
    Eigen::MatrixXd fd_jacobian;
    VectorMax12d x = pv.dof(V, E, F);
    fd::finite_jacobian(
        x,
        [&pv](const Eigen::MatrixXd& x_fd) { return pv.compute_normal(x_fd); },
        fd_jacobian);
    CHECK(fd::compare_jacobian(jacobian, fd_jacobian));
    if (!fd::compare_jacobian(jacobian, fd_jacobian)) {
        std::cout << "Jacobian:\n" << jacobian << std::endl;
        std::cout << "FD Jacobian:\n" << fd_jacobian << std::endl;
    }
}
