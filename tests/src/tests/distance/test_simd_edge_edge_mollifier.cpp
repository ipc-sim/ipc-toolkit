#include <catch2/catch_test_macros.hpp>

#include <ipc/utils/simd.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <ipc/distance/edge_edge_mollifier.hpp>

#include <Eigen/Geometry>

#include <cmath>
#include <string>
#include <vector>

using namespace ipc;

namespace {

using Batch = SimdBatch<double>;
constexpr int L = int(Batch::size);

/// @brief The threshold of two unit-length rest edges.
constexpr double EPS_X = 1e-3;

/// @brief The batch and scalar paths differ only in how the compiler contracts
/// multiply-adds. As in the other SIMD tests, the tolerance is relative to the
/// magnitude of the whole result: the mollifier derivatives scale as 1/eps_x
/// and 1/eps_x², so a per-entry absolute tolerance would be meaningless.
constexpr double TOL = 1e-9;

/// @brief Angles between the two edges placing a lane on either side of the
/// mollifier threshold: exactly parallel, mollified, just inside the
/// threshold, and inactive. The mollifier activates at sin²θ < eps_x, i.e.
/// θ < 0.0316 rad for unit rest edges.
const std::vector<double> THETAS = { 0.0, 0.02, 0.0316, 0.4 };

/// @brief Fill the lanes from `THETAS` starting at `offset`. A batch may hold
/// fewer lanes than there are cases, so the caller sweeps the offset to reach
/// every one.
std::vector<double> lane_thetas(const int offset)
{
    std::vector<double> thetas(L);
    for (int l = 0; l < L; ++l) {
        thetas[l] = THETAS[(l + offset) % THETAS.size()];
    }
    return thetas;
}

/// @brief A lane-dependent rotation, so the gradients are not mostly zeros the
/// way an axis-aligned configuration would make them.
Eigen::Matrix3d rotation(const int l)
{
    return (Eigen::AngleAxisd(
                0.7 + 0.31 * l, Eigen::Vector3d(1, 2, 3).normalized())
            * Eigen::AngleAxisd(
                0.4 * l - 0.2, Eigen::Vector3d(-2, 1, 0.5).normalized()))
        .toRotationMatrix();
}

struct Config {
    Eigen::Vector3d ea0, ea1, eb0, eb1;
};

/// @brief Two unit-length edges meeting at an angle `theta`, rigidly placed by
/// lane. Unit lengths keep the rest threshold at exactly `EPS_X`, so `theta`
/// alone decides which side of the branch the lane lands on.
Config config(const double theta, const int l)
{
    const Eigen::Matrix3d R = rotation(l);
    const Eigen::Vector3d t(0.1 * l, -0.2, 0.3);

    const Eigen::Vector3d ea0(0, 0, 0);
    const Eigen::Vector3d ea1(1, 0, 0);
    const Eigen::Vector3d eb0(0, 1, 0.3);
    const Eigen::Vector3d eb1 =
        eb0 + Eigen::Vector3d(std::cos(theta), std::sin(theta), 0);

    return { R * ea0 + t, R * ea1 + t, R * eb0 + t, R * eb1 + t };
}

/// @brief Pack the k-th coordinate of L problems into one batch.
Eigen::Vector3<Batch> pack(const std::vector<Eigen::Vector3d>& v)
{
    Eigen::Vector3<Batch> out;
    for (int k = 0; k < 3; ++k) {
        std::vector<double> tmp(L);
        for (int l = 0; l < L; ++l) {
            tmp[l] = v[l][k];
        }
        out[k] = Batch::load_unaligned(tmp.data());
    }
    return out;
}

void check_lane(
    const std::string& name,
    const int offset,
    const int l,
    const double got,
    const double want)
{
    const double scale = std::max(1.0, std::abs(want));
    CAPTURE(name, offset, l, got, want, scale);
    CHECK(std::abs(got - want) <= TOL * scale);
}

/// @brief Compare a batch-valued matrix against the scalar result per lane.
template <typename BatchMatrix, typename ScalarOf>
void check_lanes(
    const std::string& name,
    const int offset,
    const BatchMatrix& batched,
    ScalarOf&& scalar_of)
{
    for (int l = 0; l < L; ++l) {
        const auto want = scalar_of(l);
        REQUIRE(batched.rows() == want.rows());
        REQUIRE(batched.cols() == want.cols());
        // One scale for the whole result: entries of a mollifier Hessian are
        // cancellations of much larger terms.
        const double scale = std::max(1.0, want.array().abs().maxCoeff());
        for (Eigen::Index i = 0; i < want.size(); ++i) {
            CAPTURE(name, offset, l, i, scale);
            CHECK(std::abs(batched(i).get(l) - want(i)) <= TOL * scale);
        }
    }
}

} // namespace

TEST_CASE(
    "SIMD batch edge-edge mollifier matches the scalar one lane-wise",
    "[distance][mollifier][simd]")
{
    // x straddling eps_x, so one batch carries mollified and inactive lanes.
    const std::vector<double> region_xs = { 0.0, 0.1 * EPS_X, 0.9 * EPS_X,
                                            2.0 * EPS_X };

    for (int offset = 0; offset < int(region_xs.size()); ++offset) {
        std::vector<double> xs(L);
        for (int l = 0; l < L; ++l) {
            xs[l] = region_xs[(l + offset) % region_xs.size()];
        }

        const Batch x = Batch::load_unaligned(xs.data());
        const Batch eps_x(EPS_X);

        const Batch m = edge_edge_mollifier(x, eps_x);
        const Batch dm = edge_edge_mollifier_gradient(x, eps_x);
        const Batch d2m = edge_edge_mollifier_hessian(x, eps_x);
        const Batch dm_deps =
            edge_edge_mollifier_derivative_wrt_eps_x(x, eps_x);
        const Batch d2m_deps =
            edge_edge_mollifier_gradient_derivative_wrt_eps_x(x, eps_x);

        for (int l = 0; l < L; ++l) {
            check_lane(
                "mollifier", offset, l, m.get(l),
                edge_edge_mollifier(xs[l], EPS_X));
            check_lane(
                "gradient", offset, l, dm.get(l),
                edge_edge_mollifier_gradient(xs[l], EPS_X));
            check_lane(
                "hessian", offset, l, d2m.get(l),
                edge_edge_mollifier_hessian(xs[l], EPS_X));
            check_lane(
                "derivative_wrt_eps_x", offset, l, dm_deps.get(l),
                edge_edge_mollifier_derivative_wrt_eps_x(xs[l], EPS_X));
            check_lane(
                "gradient_derivative_wrt_eps_x", offset, l, d2m_deps.get(l),
                edge_edge_mollifier_gradient_derivative_wrt_eps_x(
                    xs[l], EPS_X));
        }
    }
}

TEST_CASE(
    "SIMD batch edge-edge mollifier blends the threshold branch per lane",
    "[distance][mollifier][simd]")
{
    for (int offset = 0; offset < int(THETAS.size()); ++offset) {
        const std::vector<double> thetas = lane_thetas(offset);

        std::vector<Config> configs(L);
        std::vector<Eigen::Vector3d> EA0(L), EA1(L), EB0(L), EB1(L);
        // The rest configuration is the same geometry at theta = 0: both edges
        // are unit length, so the threshold is exactly EPS_X on every lane.
        std::vector<Eigen::Vector3d> RA0(L), RA1(L), RB0(L), RB1(L);
        for (int l = 0; l < L; ++l) {
            const Config c = config(thetas[l], l);
            EA0[l] = c.ea0, EA1[l] = c.ea1, EB0[l] = c.eb0, EB1[l] = c.eb1;
            const Config r = config(0.0, l + 1);
            RA0[l] = r.ea0, RA1[l] = r.ea1, RB0[l] = r.eb0, RB1[l] = r.eb1;
        }

        const Eigen::Vector3<Batch> ea0 = pack(EA0), ea1 = pack(EA1),
                                    eb0 = pack(EB0), eb1 = pack(EB1);
        const Eigen::Vector3<Batch> ra0 = pack(RA0), ra1 = pack(RA1),
                                    rb0 = pack(RB0), rb1 = pack(RB1);
        const Batch eps_x(EPS_X);

        // The threshold itself must agree before anything built on it does.
        const Batch eps_x_batch =
            edge_edge_mollifier_threshold(ra0, ra1, rb0, rb1);
        const Batch x_batch = edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);
        const Batch m_batch = edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x);
        for (int l = 0; l < L; ++l) {
            check_lane(
                "threshold", offset, l, eps_x_batch.get(l),
                edge_edge_mollifier_threshold(RA0[l], RA1[l], RB0[l], RB1[l]));
            check_lane(
                "cross_squarednorm", offset, l, x_batch.get(l),
                edge_edge_cross_squarednorm(EA0[l], EA1[l], EB0[l], EB1[l]));
            check_lane(
                "mollifier", offset, l, m_batch.get(l),
                edge_edge_mollifier(EA0[l], EA1[l], EB0[l], EB1[l], EPS_X));
        }

        check_lanes(
            "cross_squarednorm_gradient", offset,
            edge_edge_cross_squarednorm_gradient(ea0, ea1, eb0, eb1),
            [&](int l) {
                return edge_edge_cross_squarednorm_gradient(
                    EA0[l], EA1[l], EB0[l], EB1[l]);
            });
        check_lanes(
            "cross_squarednorm_hessian", offset,
            edge_edge_cross_squarednorm_hessian(ea0, ea1, eb0, eb1),
            [&](int l) {
                return edge_edge_cross_squarednorm_hessian(
                    EA0[l], EA1[l], EB0[l], EB1[l]);
            });

        check_lanes(
            "mollifier_gradient", offset,
            edge_edge_mollifier_gradient(ea0, ea1, eb0, eb1, eps_x),
            [&](int l) {
                return edge_edge_mollifier_gradient(
                    EA0[l], EA1[l], EB0[l], EB1[l], EPS_X);
            });
        check_lanes(
            "mollifier_hessian", offset,
            edge_edge_mollifier_hessian(ea0, ea1, eb0, eb1, eps_x), [&](int l) {
                return edge_edge_mollifier_hessian(
                    EA0[l], EA1[l], EB0[l], EB1[l], EPS_X);
            });

        check_lanes(
            "threshold_gradient", offset,
            edge_edge_mollifier_threshold_gradient(ra0, ra1, rb0, rb1),
            [&](int l) {
                return edge_edge_mollifier_threshold_gradient(
                    RA0[l], RA1[l], RB0[l], RB1[l]);
            });
        check_lanes(
            "mollifier_gradient_wrt_x", offset,
            edge_edge_mollifier_gradient_wrt_x(
                ra0, ra1, rb0, rb1, ea0, ea1, eb0, eb1),
            [&](int l) {
                return edge_edge_mollifier_gradient_wrt_x(
                    RA0[l], RA1[l], RB0[l], RB1[l], EA0[l], EA1[l], EB0[l],
                    EB1[l]);
            });
        check_lanes(
            "mollifier_gradient_jacobian_wrt_x", offset,
            edge_edge_mollifier_gradient_jacobian_wrt_x(
                ra0, ra1, rb0, rb1, ea0, ea1, eb0, eb1),
            [&](int l) {
                return edge_edge_mollifier_gradient_jacobian_wrt_x(
                    RA0[l], RA1[l], RB0[l], RB1[l], EA0[l], EA1[l], EB0[l],
                    EB1[l]);
            });
    }
}

#endif
