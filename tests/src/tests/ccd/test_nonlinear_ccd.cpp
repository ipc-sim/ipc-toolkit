#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/ccd/nonlinear_ccd.hpp>
#include <ipc/distance/point_line.hpp>

#include <igl/PI.h>

using namespace ipc;

class RotationalTrajectory : public NonlinearTrajectory {
public:
    RotationalTrajectory(
        const VectorMax3d& point,
        const VectorMax3d& center,
        const double z_angular_velocity)
        : center(center)
        , point(point)
        , z_angular_velocity(z_angular_velocity)
    {
        assert(point.size() == point.size());
    }

    VectorMax3d operator()(const double t) const override
    {
        return position(t);
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

private:
    template <typename T> VectorMax3<T> position(const T& t) const
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

class StaticTrajectory : public NonlinearTrajectory {
public:
    StaticTrajectory(const VectorMax3d& point) : point(point) { }

    VectorMax3d operator()(const double t) const override { return point; }

    double
    max_distance_from_linear(const double t0, const double t1) const override
    {
        return 0;
    }

protected:
    VectorMax3d point;
};

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
