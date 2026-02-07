#include <Eigen/Core>
#include <catch2/catch_all.hpp>

#include <ipc/dynamics/rigid/pose.hpp>

#include <Eigen/Geometry>
#include <catch2/generators/catch_generators_range.hpp>
#include <igl/PI.h>
#include <finitediff.hpp>

#include <iostream>

using namespace ipc;
using namespace ipc::rigid;

TEST_CASE("so(3) -> SO(3)", "[rigid][pose]")
{
    double angle;
    Eigen::Vector3d axis;

    SECTION("zero")
    {
        angle = 0;
        axis = Eigen::Vector3d::Random();
    }
    SECTION("small angle")
    {
        angle = GENERATE(take(100, random(0.0, 1e-6)));
        axis = Eigen::Vector3d::Random();
    }
    SECTION("random")
    {
        angle = GENERATE(take(100, random(0.0, 2 * igl::PI)));
        axis = Eigen::Vector3d::Random();
    }
    axis.normalize();

    const Eigen::Matrix3d R_expected =
        Eigen::AngleAxisd(angle, axis).toRotationMatrix();

    const Eigen::Vector3d theta = angle * axis;
    const Eigen::Matrix3d R_actual = rotation_vector_to_matrix(theta);

    CHECK((R_actual - R_expected).norm() == Catch::Approx(0).margin(1e-12));
}

TEST_CASE("Benchmark so(3) -> SO(3)", "[!benchmark][rigid][pose]")
{
    Eigen::Vector3d theta = Eigen::Vector3d::Random();

    BENCHMARK("Compute R") { return rotation_vector_to_matrix(theta); };

    BENCHMARK("Compute ∇R")
    {
        return rotation_vector_to_matrix_jacobian(theta);
    };

    BENCHMARK("Compute ∇²R")
    {
        return rotation_vector_to_matrix_hessian(theta);
    };
}

#if false
TEST_CASE("Interval so(3) -> SO(3)", "[!benchmark][physics][pose]")
{
    using namespace ipc::rigid;
    double angle;
    Eigen::Vector3d axis;

    SECTION("zero")
    {
        angle = 0;
        axis = Eigen::Vector3d::Random();
    }
    SECTION("random")
    {
        angle = GENERATE(take(1, random(0.0, 2 * igl::PI)));
        axis = Eigen::Vector3d::Random();
    }
    axis.normalize();

    Pose<double> p = Pose<double>::Zero(3);
    p.rotation = angle * axis;
    BENCHMARK("Double so(3) -> SO(3)")
    {
        Eigen::Matrix3d R = p.construct_rotation_matrix();
    };
    Pose<ipc::rigid::Interval> pI = p.cast<ipc::rigid::Interval>();
    BENCHMARK("Interval so(3) -> SO(3)")
    {
        Matrix3I R = pI.construct_rotation_matrix();
    };
}
#endif

TEST_CASE("so(3) -> SO(3) derivatives", "[rigid][pose]")
{
    double angle;
    Eigen::Vector3d axis;

    SECTION("zero")
    {
        angle = 0;
        axis = Eigen::Vector3d::Random();
    }
    SECTION("small angle")
    {
        angle = GENERATE(take(100, random(0.0, 1e-6)));
        axis = Eigen::Vector3d::Random();
    }
    SECTION("random")
    {
        angle = GENERATE(take(100, random(0.0, 2 * igl::PI)));
        axis = Eigen::Vector3d::Random();
    }
    axis.normalize();

    Eigen::Vector3d theta = angle * axis;

    {
        Eigen::Matrix<double, 9, 3> J_analytic =
            rotation_vector_to_matrix_jacobian(theta);

        auto f = [&](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
            return rotation_vector_to_matrix(x);
        };

        Eigen::MatrixXd J_numerical;
        fd::finite_jacobian_tensor<3>(theta, f, J_numerical);

        CHECK(fd::compare_jacobian(J_analytic, J_numerical));
        if (!fd::compare_jacobian(J_analytic, J_numerical)) {
            std::cout << "J_analytic:\n" << J_analytic << "\n\n";
            std::cout << "J_numerical:\n" << J_numerical << "\n\n";
        }
    }

    {
        Eigen::Matrix<double, 9, 9> H_analytic =
            rotation_vector_to_matrix_hessian(theta);

        auto f = [&](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
            return rotation_vector_to_matrix_jacobian(x);
        };

        Eigen::MatrixXd H_numerical;
        fd::finite_jacobian_tensor<4>(theta, f, H_numerical);

        CHECK(fd::compare_jacobian(H_analytic, H_numerical));
        if (!fd::compare_jacobian(H_analytic, H_numerical)) {
            std::cout << "H_analytic:\n" << H_analytic << "\n\n";
            std::cout << "H_numerical:\n" << H_numerical << "\n\n";
        }
    }
}

TEST_CASE("SO(3) -> so(3)", "[rigid][pose]")
{
    double angle = GENERATE(
        0.0, 1e-8, 1e-6, igl::PI / 4, igl::PI / 2, igl::PI * 0.75,
        igl::PI - 1e-4, igl::PI - 1e-6, igl::PI - 1e-8, igl::PI, igl::PI + 1e-8,
        igl::PI + 1e-6, igl::PI + 1e-4, igl::PI * 1.5, igl::PI * 2);
    Eigen::Vector3d axis = Eigen::Vector3d::Random().normalized();

    Eigen::Vector3d theta_expected = angle * axis;

    Eigen::Matrix3d R = rotation_vector_to_matrix(theta_expected);

    Eigen::Vector3d theta = rotation_matrix_to_vector(R);

    Eigen::AngleAxisd angle_axis(R);
    theta_expected = angle_axis.angle() * angle_axis.axis();

    CHECK(
        (theta.isApprox(theta_expected, 1e-4)
         || theta.isApprox(angle * axis, 1e-4)));

    // std::cout << "input axis: " << axis.transpose() << "\n"
    //           << "input angle: " << angle << " (" << angle - igl::PI << ")"
    //           << "\n"
    //           << "expected axis: " << angle_axis.axis().transpose() << "\n"
    //           << "expected angle: " << angle_axis.angle() << "\n"
    //           << "axis: " << theta.normalized().transpose() << "\n"
    //           << "angle: " << theta.norm() << "\n"
    //           << "\n"
    //           << std::endl;
}

TEST_CASE("Benchmark SO(3) -> so(3)", "[!benchmark][rigid][pose]")
{
    Eigen::Matrix3d R = Eigen::Matrix3d::Random();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(
        R, Eigen::ComputeFullU | Eigen::ComputeFullV);
    R = svd.matrixU() * svd.matrixV().transpose(); // closest rotation matrix

    BENCHMARK("Mine") { return rotation_matrix_to_vector(R); };

    BENCHMARK("Eigen")
    {
        Eigen::AngleAxisd angle_axis(R);
        return angle_axis.angle() * angle_axis.axis();
    };
}

TEST_CASE("Pose transform vertices Jacobian", "[rigid][pose][jacobian]")
{
    GENERATE(range(0, 10));
    const int DIM = GENERATE(2, 3);

    Pose p;
    p.position = VectorMax3d::Random(DIM);
    p.rotation = VectorMax3d::Random(DIM == 2 ? 1 : 3);

    constexpr int NUM_VERTICES = 2;
    Eigen::MatrixXd V = Eigen::MatrixXd::Random(NUM_VERTICES, DIM);

    Eigen::MatrixXd dV_dx_analytic = p.transform_vertices_jacobian(V);

    auto f = [&](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
        return Pose(x).transform_vertices(V);
    };

    Eigen::VectorXd x(DIM == 2 ? 3 : 6);
    x << p.position, p.rotation;

    Eigen::MatrixXd dV_dx_numerical;
    fd::finite_jacobian_tensor<3>(x, f, dV_dx_numerical);

    CHECK(fd::compare_jacobian(dV_dx_analytic, dV_dx_numerical));
    if (!fd::compare_jacobian(dV_dx_analytic, dV_dx_numerical)) {
        std::cout << "dV_dx_analytic:\n" << dV_dx_analytic << "\n\n";
        std::cout << "dV_dx_numerical:\n" << dV_dx_numerical << "\n\n";
    }
}

TEST_CASE("Pose transform vertices Hessian", "[rigid][pose][hessian]")
{
    // GENERATE(range(0, 10));
    const int DIM = GENERATE(2, 3);

    Pose p;
    p.position = VectorMax3d::Random(DIM);
    p.rotation = VectorMax3d::Random(DIM == 2 ? 1 : 3);

    constexpr int NUM_VERTICES = 2;
    Eigen::MatrixXd V = Eigen::MatrixXd::Random(NUM_VERTICES, DIM);

    Eigen::MatrixXd d2V_dx2_analytic = p.transform_vertices_hessian(V);

    auto f = [&](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
        return Pose(x).transform_vertices_jacobian(V);
    };

    Eigen::VectorXd x(DIM == 2 ? 3 : 6);
    x << p.position, p.rotation;

    Eigen::MatrixXd d2V_dx2_numerical;
    fd::finite_jacobian_tensor<4>(x, f, d2V_dx2_numerical);

    CHECK(fd::compare_jacobian(d2V_dx2_analytic, d2V_dx2_numerical));
    if (!fd::compare_jacobian(d2V_dx2_analytic, d2V_dx2_numerical)) {
        std::cout << "dV_dx_analytic:\n" << d2V_dx2_analytic << "\n\n";
        std::cout << "dV_dx_numerical:\n" << d2V_dx2_numerical << "\n\n";
    }
}