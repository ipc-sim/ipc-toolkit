#include <catch2/catch_all.hpp>

#include <ipc/dynamics/rigid/pose.hpp>

#include <Eigen/Geometry>
#include <igl/PI.h>

using namespace ipc::rigid;

TEST_CASE("SE(3) ↦ SO(3)", "[rigid][pose]")
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

TEST_CASE("∇²(SE(3) ↦ SO(3))", "[!benchmark][rigid][pose]")
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
TEST_CASE("Interval SE(3) ↦ SO(3)", "[!benchmark][physics][pose]")
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
    BENCHMARK("Double SE(3) ↦ SO(3)")
    {
        Eigen::Matrix3d R = p.construct_rotation_matrix();
    };
    Pose<ipc::rigid::Interval> pI = p.cast<ipc::rigid::Interval>();
    BENCHMARK("Interval SE(3) ↦ SO(3)")
    {
        Matrix3I R = pI.construct_rotation_matrix();
    };
}
#endif
