#include <Eigen/Core>
#include <catch2/catch_all.hpp>

#include <ipc/ccd/nonlinear_ccd.hpp>
#include <ipc/dynamics/rigid/rigid_trajectory.hpp>
#ifdef IPC_TOOLKIT_WITH_FILIB
#include <ipc/math/interval.hpp>
#include <ipc/math/sinc.hpp>
#endif

#include <igl/PI.h>

using namespace ipc;
using namespace ipc::rigid;

namespace {

/// Densely sample max_t ‖x(t) − lerp(x(t0), x(t1))(t)‖ over [t0, t1].
double sampled_max_deviation(
    const RigidTrajectory& traj,
    const double t0,
    const double t1,
    const int n = 1000)
{
    const VectorMax3d x_t0 = traj(t0);
    const VectorMax3d x_t1 = traj(t1);
    double max_d = 0;
    for (int i = 0; i <= n; i++) {
        const double s = i / double(n);
        const double t = t0 + s * (t1 - t0);
        max_d = std::max(max_d, (traj(t) - (x_t0 + s * (x_t1 - x_t0))).norm());
    }
    return max_d;
}

Pose random_pose(const int dim, const double max_angle)
{
    VectorMax3d position = VectorMax3d::Random(dim);
    VectorMax3d rotation;
    if (dim == 2) {
        rotation = max_angle * VectorMax3d::Random(1);
    } else {
        rotation = VectorMax3d::Random(3).normalized()
            * (max_angle * 0.5 * (Eigen::Vector2d::Random()(0) + 1));
    }
    return Pose(position, rotation);
}

} // namespace

TEST_CASE("Rigid trajectory position", "[rigid][trajectory]")
{
    const int dim = GENERATE(2, 3);
    for (int trial = 0; trial < 10; trial++) {
        const Pose pose_t0 = random_pose(dim, 2 * igl::PI);
        const Pose pose_t1 = random_pose(dim, 2 * igl::PI);
        const VectorMax3d rest_position = VectorMax3d::Random(dim);

        const RigidTrajectory traj(rest_position, pose_t0, pose_t1);

        for (int i = 0; i <= 100; i++) {
            const double t = i / 100.0;
            const Pose pose_t(
                (1 - t) * pose_t0.position + t * pose_t1.position,
                (1 - t) * pose_t0.rotation + t * pose_t1.rotation);
            const Eigen::MatrixXd expected =
                pose_t.transform_vertices(rest_position.transpose());
            CHECK(
                (traj(t) - expected.row(0).transpose()).norm()
                == Catch::Approx(0).margin(1e-14));
        }

        // Endpoints are exact
        CHECK(
            (traj(0)
             - (pose_t0.rotation_matrix() * rest_position + pose_t0.position))
                .norm()
            == Catch::Approx(0).margin(1e-15));
        CHECK(
            (traj(1)
             - (pose_t1.rotation_matrix() * rest_position + pose_t1.position))
                .norm()
            == Catch::Approx(0).margin(1e-15));
    }
}

TEST_CASE(
    "Rigid trajectory max_distance_from_linear is conservative",
    "[rigid][trajectory]")
{
    const int dim = GENERATE(2, 3);
    double max_angle;
    SECTION("small rotation") { max_angle = 0.1; }
    SECTION("moderate rotation") { max_angle = igl::PI / 2; }
    SECTION("large rotation") { max_angle = 2 * igl::PI; }

    for (int trial = 0; trial < 10; trial++) {
        const Pose pose_t0 = random_pose(dim, max_angle);
        const Pose pose_t1 = random_pose(dim, max_angle);
        const VectorMax3d rest_position = VectorMax3d::Random(dim);

        const RigidTrajectory traj(rest_position, pose_t0, pose_t1);

        // Check on the full interval and random sub-intervals
        std::vector<std::pair<double, double>> intervals = { { 0, 1 } };
        for (int i = 0; i < 4; i++) {
            const Eigen::Vector2d r =
                (Eigen::Vector2d::Random() + Eigen::Vector2d::Ones()) / 2;
            intervals.emplace_back(
                std::min(r(0), r(1)), std::max(r(0), r(1)));
        }

        for (const auto& [t0, t1] : intervals) {
            if (t1 - t0 < 1e-8) {
                continue;
            }
            const double bound = traj.max_distance_from_linear(t0, t1);
            const double sampled = sampled_max_deviation(traj, t0, t1);
            CHECK(sampled <= bound * (1 + 1e-10) + 1e-12);
        }
    }
}

TEST_CASE(
    "Rigid trajectory pure translation has zero deviation",
    "[rigid][trajectory]")
{
    const int dim = GENERATE(2, 3);
    const VectorMax3d rotation =
        dim == 2 ? VectorMax3d::Random(1) : VectorMax3d::Random(3);
    const Pose pose_t0(VectorMax3d::Random(dim), rotation);
    const Pose pose_t1(VectorMax3d::Random(dim), rotation); // same rotation

    const RigidTrajectory traj(VectorMax3d::Random(dim), pose_t0, pose_t1);

    CHECK(traj.max_distance_from_linear(0, 1) == 0.0);
    CHECK(sampled_max_deviation(traj, 0, 1) == Catch::Approx(0).margin(1e-14));
}

TEST_CASE("Rigid trajectory bound tightness", "[rigid][trajectory]")
{
    // Fixed-axis rotation with the point perpendicular to the axis: the true
    // deviation is the sagitta r(1 − cos(δ/2)); the bound is r δ²/8. Check the
    // bound is not grossly over-conservative in this regime.
    const double delta = GENERATE(0.1, 0.5, 1.0, igl::PI / 2);

    const Pose pose_t0(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
    const Pose pose_t1(
        Eigen::Vector3d::Zero(), Eigen::Vector3d(0, 0, delta));
    const RigidTrajectory traj(Eigen::Vector3d(1, 0, 0), pose_t0, pose_t1);

    const double bound = traj.max_distance_from_linear(0, 1);
    const double sampled = sampled_max_deviation(traj, 0, 1);

    CHECK(sampled <= bound);
    CHECK(bound <= 2 * sampled); // tightness (ratio → 1 as δ → 0)
}

// ============================================================================
// Interval-arithmetic cross-check (tests only; production code is IA-free)
#ifdef IPC_TOOLKIT_WITH_FILIB

namespace {

class IntervalRigidTrajectory : public RigidTrajectory,
                                public IntervalNonlinearTrajectory {
public:
    IntervalRigidTrajectory(
        Eigen::ConstRef<VectorMax3d> rest_position,
        const Pose& pose_t0,
        const Pose& pose_t1)
        : RigidTrajectory(rest_position, pose_t0, pose_t1)
        , m_rest_position(rest_position)
        , m_pose_t0(pose_t0)
        , m_pose_t1(pose_t1)
    {
    }

    using RigidTrajectory::operator();

    VectorMax3I operator()(const filib::Interval& t) const override
    {
        const filib::Interval one_minus_t = filib::Interval(1.0) - t;

        const VectorMax3I p = one_minus_t
                * m_pose_t0.position.template cast<filib::Interval>()
            + t * m_pose_t1.position.template cast<filib::Interval>();
        const VectorMax3I theta = one_minus_t
                * m_pose_t0.rotation.template cast<filib::Interval>()
            + t * m_pose_t1.rotation.template cast<filib::Interval>();

        const VectorMax3I x_rest =
            m_rest_position.template cast<filib::Interval>();

        if (theta.size() == 1) {
            // 2D rotation
            MatrixMax3I R(2, 2);
            R(0, 0) = cos(theta(0));
            R(0, 1) = -sin(theta(0));
            R(1, 0) = sin(theta(0));
            R(1, 1) = cos(theta(0));
            return R * x_rest + p;
        }

        // 3D: Rodrigues formula R = I + sinc(‖θ‖) K + ½ sinc²(‖θ‖/2) K²
        const filib::Interval angle = norm(theta);
        Matrix3I K;
        K << filib::Interval(0.0), -theta(2), theta(1), //
            theta(2), filib::Interval(0.0), -theta(0),  //
            -theta(1), theta(0), filib::Interval(0.0);
        const filib::Interval s = sinc(angle);
        const filib::Interval s_half = sinc(angle * filib::Interval(0.5));
        Matrix3I R = s * K
            + (filib::Interval(0.5) * s_half * s_half) * (K * K).eval();
        R.diagonal().array() += filib::Interval(1.0);
        return R * x_rest + p;
    }

    // Expose the interval-based estimate for comparison
    double interval_max_distance_from_linear(double t0, double t1) const
    {
        return IntervalNonlinearTrajectory::max_distance_from_linear(t0, t1);
    }

    double
    max_distance_from_linear(const double t0, const double t1) const override
    {
        return RigidTrajectory::max_distance_from_linear(t0, t1);
    }

private:
    VectorMax3d m_rest_position;
    Pose m_pose_t0, m_pose_t1;
};

} // namespace

TEST_CASE(
    "Rigid trajectory interval cross-check", "[rigid][trajectory][filib]")
{
    const int dim = GENERATE(2, 3);
    const double max_angle = GENERATE(0.1, 1.0, igl::PI);

    for (int trial = 0; trial < 5; trial++) {
        const Pose pose_t0 = random_pose(dim, max_angle);
        const Pose pose_t1 = random_pose(dim, max_angle);
        const VectorMax3d rest_position = VectorMax3d::Random(dim);

        const IntervalRigidTrajectory traj(rest_position, pose_t0, pose_t1);

        // Interval evaluation contains the true position
        for (int i = 0; i < 10; i++) {
            const double t = i / 9.0;
            const VectorMax3d x = traj(t);
            const VectorMax3I x_I = traj(filib::Interval(t));
            for (int d = 0; d < dim; d++) {
                CHECK(x(d) >= x_I(d).INF - 1e-12);
                CHECK(x(d) <= x_I(d).SUP + 1e-12);
            }
        }

        // Both the analytic bound and the interval estimate dominate the
        // sampled truth
        const double sampled = sampled_max_deviation(traj, 0, 1);
        const double analytic = traj.max_distance_from_linear(0, 1);
        const double interval = traj.interval_max_distance_from_linear(0, 1);
        CHECK(sampled <= analytic * (1 + 1e-10) + 1e-12);
        CHECK(sampled <= interval * (1 + 1e-10) + 1e-12);
    }
}

#endif

// ============================================================================
// CCD with rigid trajectories (mirrors test_nonlinear_ccd.cpp scenarios)

TEST_CASE("Rigid trajectory point-edge CCD", "[rigid][trajectory][ccd]")
{
    // Static point at (0, 0.5); edge from (-1,0) to (1,0) rotating by π.
    const Pose identity_2d = Pose::Identity(2);
    const Pose rotated_2d(
        Eigen::Vector2d::Zero(), VectorMax3d::Constant(1, igl::PI));

    const RigidTrajectory p(
        Eigen::Vector2d(0, 0.5), identity_2d, identity_2d);
    const RigidTrajectory e0(Eigen::Vector2d(-1, 0), identity_2d, rotated_2d);
    const RigidTrajectory e1(Eigen::Vector2d(1, 0), identity_2d, rotated_2d);

    NonlinearCCD ccd;
    ccd.conservative_rescaling = 0.9;

    double toi;
    const bool collision = ccd.point_edge_ccd(p, e0, e1, toi);

    CHECK(collision);
    CHECK((0.49 <= toi && toi <= 0.5)); // conservative estimate
}

TEST_CASE("Rigid trajectory edge-edge CCD", "[rigid][trajectory][ccd]")
{
    // Edge from (-1,0,0) to (1,0,0) rotating 2π about z; static edge at
    // y = 0.5. First contact at 30° of rotation.
    const Pose identity_3d = Pose::Identity(3);
    const Pose rotated_3d(
        Eigen::Vector3d::Zero(), Eigen::Vector3d(0, 0, 2 * igl::PI));

    const RigidTrajectory ea0(
        Eigen::Vector3d(-1, 0, 0), identity_3d, rotated_3d);
    const RigidTrajectory ea1(
        Eigen::Vector3d(1, 0, 0), identity_3d, rotated_3d);
    const RigidTrajectory eb0(
        Eigen::Vector3d(-1, 0.5, 0), identity_3d, identity_3d);
    const RigidTrajectory eb1(
        Eigen::Vector3d(1, 0.5, 0), identity_3d, identity_3d);

    double toi;
    const bool collision = NonlinearCCD().edge_edge_ccd(ea0, ea1, eb0, eb1, toi);

    CHECK(collision);
    CHECK(toi <= 30 / 360.0);
    CHECK(toi == Catch::Approx(30 / 360.0).margin(1e-2));
}

TEST_CASE("Rigid trajectory point-triangle CCD", "[rigid][trajectory][ccd]")
{
    const double x = GENERATE(-0.1, 0.0, 0.1);

    const Pose identity_3d = Pose::Identity(3);
    const Pose rotated_3d(
        Eigen::Vector3d::Zero(), Eigen::Vector3d(0, 0, igl::PI));

    const RigidTrajectory t0(
        Eigen::Vector3d(1, 0, 0), identity_3d, rotated_3d);
    const RigidTrajectory t1(
        Eigen::Vector3d(x, 0, 1), identity_3d, rotated_3d);
    const RigidTrajectory t2(
        Eigen::Vector3d(x, 0, -1), identity_3d, rotated_3d);
    const RigidTrajectory p(
        Eigen::Vector3d(0, 0.5, 0), identity_3d, identity_3d);

    NonlinearCCD ccd;
    ccd.conservative_rescaling = 0.9;

    double toi;
    const bool collision = ccd.point_triangle_ccd(p, t0, t1, t2, toi);

    CHECK(collision);
    CHECK(toi <= 0.5);
    CHECK(toi == Catch::Approx(0.5).margin(1e-2));
}
