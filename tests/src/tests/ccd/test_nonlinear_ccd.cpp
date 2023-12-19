#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/ccd/nonlinear_ccd.hpp>
#include <ipc/distance/point_line.hpp>

#include <igl/PI.h>

using namespace ipc;

class RotationalTrajectory : virtual public NonlinearTrajectory {
public:
    RotationalTrajectory(
        const VectorMax3d& _point,
        const VectorMax3d& _center,
        const double _z_angular_velocity)
        : center(_center)
        , point(_point)
        , z_angular_velocity(_z_angular_velocity)
    {
        assert(point.size() == point.size());
    }

    VectorMax3d operator()(const double t) const override
    {
        return position(center, point, z_angular_velocity, t);
    }

    double
    max_distance_from_linear(const double t0, const double t1) const override
    {
        if (z_angular_velocity * (t1 - t0) >= 2 * igl::PI) {
            // This is the most conservative estimate
            return 2 * (point - center).norm(); // 2 * radius
        }

        const VectorMax3d p_t0 = (*this)(t0);
        const VectorMax3d p_t1 = (*this)(t1);
        return ((*this)((t0 + t1) / 2) - ((p_t1 - p_t0) * 0.5 + p_t0)).norm();
    }

protected:
    VectorMax3d center, point;
    double z_angular_velocity;

    template <typename T>
    static VectorMax3<T> position(
        const VectorMax3d& center,
        const VectorMax3d& point,
        const double z_angular_velocity,
        const T& t)
    {
        MatrixMax3<T> R(point.size(), point.size());

        R(0, 0) = cos(z_angular_velocity * t);
        R(0, 1) = -sin(z_angular_velocity * t);
        R(1, 0) = sin(z_angular_velocity * t);
        R(1, 1) = cos(z_angular_velocity * t);
        if (point.size() == 3) {
            R(0, 2) = R(1, 2) = R(2, 0) = R(2, 1) = T(0);
            R(2, 2) = T(1);
        }

        return R * (point - center) + center;
    }
};

class IntervalRotationalTrajectory : public RotationalTrajectory,
                                     public IntervalNonlinearTrajectory {
public:
    IntervalRotationalTrajectory(
        const VectorMax3d& _point,
        const VectorMax3d& _center,
        const double _z_angular_velocity)
        : RotationalTrajectory(_point, _center, _z_angular_velocity)
    {
    }

    using RotationalTrajectory::operator();

    VectorMax3I operator()(const filib::Interval& t) const override
    {
        return position(center, point, z_angular_velocity, t);
    }

    double
    max_distance_from_linear(const double t0, const double t1) const override
    {
        return IntervalNonlinearTrajectory::max_distance_from_linear(t0, t1);
    }
};

class StaticTrajectory : public NonlinearTrajectory {
public:
    StaticTrajectory(const VectorMax3d& _point) : point(_point) { }

    VectorMax3d operator()(const double t) const override { return point; }

    double
    max_distance_from_linear(const double t0, const double t1) const override
    {
        return 0;
    }

protected:
    VectorMax3d point;
};

// BEGIN_RIGID_2D_TRAJECTORY
class Rigid2DTrajectory : virtual public ipc::NonlinearTrajectory {
public:
    Rigid2DTrajectory(
        const Eigen::Vector2d& _position,
        const Eigen::Vector2d& _translation,
        const Eigen::Vector2d& _delta_translation,
        const double _rotation,
        const double _delta_rotation)
        : position(_position)
        , translation(_translation)
        , delta_translation(_delta_translation)
        , rotation(_rotation)
        , delta_rotation(_delta_rotation)
    {
    }

    VectorMax3d operator()(const double t) const override
    {
        const Eigen::Matrix2d R =
            Eigen::Rotation2D<double>(rotation + t * delta_rotation)
                .toRotationMatrix();

        return R * position + translation + t * delta_translation;
    }

    double
    max_distance_from_linear(const double t0, const double t1) const override
    {
        if (delta_rotation * (t1 - t0) >= 2 * igl::PI) {
            // This is the most conservative estimate
            return 2 * position.norm(); // 2 * radius
        }

        const VectorMax3d p_t0 = (*this)(t0);
        const VectorMax3d p_t1 = (*this)(t1);
        return ((*this)((t0 + t1) / 2) - ((p_t1 - p_t0) * 0.5 + p_t0)).norm();
    }

protected:
    Eigen::Vector2d position;
    Eigen::Vector2d translation;
    Eigen::Vector2d delta_translation;
    double rotation;
    double delta_rotation;
};
// END_RIGID_2D_TRAJECTORY

TEST_CASE("Nonlinear Point-Point CCD", "[ccd][nonlinear][point-point]")
{
    const StaticTrajectory p0(Eigen::Vector2d(0, 1));

    std::unique_ptr<RotationalTrajectory> p1;
    SECTION("Analytic max")
    {
        p1 = std::make_unique<RotationalTrajectory>(
            Eigen::Vector2d(1, 0), Eigen::Vector2d::Zero(), igl::PI);
    }
    SECTION("Interval max")
    {
        p1 = std::make_unique<IntervalRotationalTrajectory>(
            Eigen::Vector2d(1, 0), Eigen::Vector2d::Zero(), igl::PI);
    }

    double toi;
    bool collision = point_point_nonlinear_ccd(
        p0, *p1, toi, /*tmax=*/1.0, /*min_distance=*/0, DEFAULT_CCD_TOLERANCE,
        DEFAULT_CCD_MAX_ITERATIONS, /*conservative_rescaling=*/0.9);

    CHECK(collision);
    CHECK(toi <= 0.5);
    CHECK(toi == Catch::Approx(0.5).margin(1e-2));
}

TEST_CASE("Nonlinear Point-Edge CCD", "[ccd][nonlinear][point-edge]")
{
    const StaticTrajectory p(Eigen::Vector2d(0, 0.5));
    const RotationalTrajectory e0(
        Eigen::Vector2d(-1, 0), Eigen::Vector2d::Zero(), 2 * igl::PI);
    const RotationalTrajectory e1(
        Eigen::Vector2d(1, 0), Eigen::Vector2d::Zero(), 2 * igl::PI);

    double toi;
    bool collision = point_edge_nonlinear_ccd(
        p, e0, e1, toi, /*tmax=*/1.0, /*min_distance=*/0, DEFAULT_CCD_TOLERANCE,
        DEFAULT_CCD_MAX_ITERATIONS, /*conservative_rescaling=*/0.9);

    CHECK(collision);
    CHECK(toi <= 0.25);
    CHECK(toi == Catch::Approx(0.25).margin(1e-2));
}

TEST_CASE("Rigid 2D Trajectory", "[ccd][nonlinear][point-edge]")
{
    // BEGIN_TEST_RIGID_2D_TRAJECTORY
    // Static point
    const Rigid2DTrajectory p(
        Eigen::Vector2d(0, 0.5), Eigen::Vector2d::Zero(),
        Eigen::Vector2d::Zero(), 0, 0);
    // Rotating edge
    const Rigid2DTrajectory e0(
        Eigen::Vector2d(-1, 0), Eigen::Vector2d::Zero(),
        Eigen::Vector2d::Zero(), 0, igl::PI);
    const Rigid2DTrajectory e1(
        Eigen::Vector2d(+1, 0), Eigen::Vector2d::Zero(),
        Eigen::Vector2d::Zero(), 0, igl::PI);

    double toi;
    bool collision = ipc::point_edge_nonlinear_ccd(
        p, e0, e1, toi, /*tmax=*/1.0, /*min_distance=*/0, DEFAULT_CCD_TOLERANCE,
        DEFAULT_CCD_MAX_ITERATIONS,
        // increase the conservative_rescaling from 0.8 to 0.9 to get a more
        // accurate estimate
        /*conservative_rescaling=*/0.9);

    CHECK(collision);
    CHECK((0.49 <= toi && toi <= 0.5)); // conservative estimate
    // END_TEST_RIGID_2D_TRAJECTORY
}

TEST_CASE("Nonlinear Edge-Edge CCD", "[ccd][nonlinear][edge-edge]")
{
    const RotationalTrajectory ea0(
        Eigen::Vector3d(-1, 0, 0), Eigen::Vector3d::Zero(), 2 * igl::PI);
    const RotationalTrajectory ea1(
        Eigen::Vector3d(1, 0, 0), Eigen::Vector3d::Zero(), 2 * igl::PI);
    const StaticTrajectory eb0(Eigen::Vector3d(-1, 0.5, 0));
    const StaticTrajectory eb1(Eigen::Vector3d(1, 0.5, 0));

    double toi;
    bool collision = edge_edge_nonlinear_ccd(ea0, ea1, eb0, eb1, toi);

    CHECK(collision);
    CHECK(toi <= 30 / 360.0);
    CHECK(toi == Catch::Approx(30 / 360.0).margin(1e-2));
}

TEST_CASE("Nonlinear Point-Triangle CCD", "[ccd][nonlinear][point-triangle]")
{
    const double x = GENERATE(-0.1, 0, 0.1);

    const RotationalTrajectory t0(
        Eigen::Vector3d(1, 0, 0), Eigen::Vector3d::Zero(), igl::PI);
    const RotationalTrajectory t1(
        Eigen::Vector3d(x, 0, 1), Eigen::Vector3d::Zero(), igl::PI);
    const RotationalTrajectory t2(
        Eigen::Vector3d(x, 0, -1), Eigen::Vector3d::Zero(), igl::PI);
    const StaticTrajectory p(Eigen::Vector3d(0, 0.5, 0));

    double toi;
    bool collision = point_triangle_nonlinear_ccd(
        p, t0, t1, t2, toi, /*tmax=*/1.0, /*min_distance=*/0,
        DEFAULT_CCD_TOLERANCE, DEFAULT_CCD_MAX_ITERATIONS,
        /*conservative_rescaling=*/0.9);

    CHECK(collision);
    CHECK(toi <= 0.5);
    CHECK(toi == Catch::Approx(0.5).margin(1e-2));
}