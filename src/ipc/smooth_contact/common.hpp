#pragma once
#include <ipc/utils/logger.hpp>

#include <cmath>

namespace ipc {

static constexpr int N_VERT_NEIGHBORS_2D = 3;
static constexpr int N_EDGE_NEIGHBORS_2D = 2;
static constexpr int MAX_VERT_2D =
    2 * std::max(N_VERT_NEIGHBORS_2D, N_EDGE_NEIGHBORS_2D);
static constexpr int N_VERT_NEIGHBORS_3D = 20; // increase me if needed
static constexpr int N_EDGE_NEIGHBORS_3D = 4;
static constexpr int N_FACE_NEIGHBORS_3D = 3;
static constexpr int MAX_VERT_3D = N_VERT_NEIGHBORS_3D * 2;

template <int DIM> class MaxVertices;
template <> class MaxVertices<2> {
public:
    static constexpr int value = MAX_VERT_2D; // NOLINT
};
template <> class MaxVertices<3> {
public:
    static constexpr int value = MAX_VERT_3D; // NOLINT
};

struct SmoothContactParameters {
    SmoothContactParameters() = default;

    SmoothContactParameters(
        const double _dhat,
        const double _alpha_t,
        const double _beta_t,
        const int _r)
        : SmoothContactParameters(_dhat, _alpha_t, _beta_t, 0, 0.1, _r)
    {
    }

    SmoothContactParameters(
        const double _dhat,
        const double _alpha_t,
        const double _beta_t,
        const double _alpha_n,
        const double _beta_n,
        const int _r)
        : dhat(_dhat)
        , alpha_t(_alpha_t)
        , beta_t(_beta_t)
        , alpha_n(_alpha_n)
        , beta_n(_beta_n)
        , r(_r)
    {
        if (r <= 0) {
            logger().error("Parameter 'r' must be greater than 0! r: {}", r);
        }
        if (dhat <= 0) {
            logger().error(
                "Parameter 'dhat' must be greater than 0! dhat: {}", dhat);
        }
        if (abs(alpha_t) > 1) {
            logger().error(
                "Parameter 'alpha_t' must be in [-1, 1]! alpha_t: {}", alpha_t);
        }
        if (abs(alpha_n) > 1) {
            logger().error(
                "Parameter 'alpha_n' must be in [-1, 1]! alpha_n: {}", alpha_n);
        }
        if (abs(beta_t) > 1) {
            logger().error(
                "Parameter 'beta_t' must be in [-1, 1]! beta_t: {}", beta_t);
        }
        if (abs(beta_n) > 1) {
            logger().error(
                "Parameter 'beta_n' must be in [-1, 1]! beta_n: {}", beta_n);
        }
        if (beta_t + alpha_t <= 1e-6) {
            logger().error(
                "Sum of 'beta_t' and 'alpha_t' must be greater than 1e-6! beta_t: {}, alpha_t: {}",
                beta_t, alpha_t);
        }
        if (beta_n + alpha_n <= 1e-6) {
            logger().error(
                "Sum of 'beta_n' and 'alpha_n' must be greater than 1e-6! beta_n: {}, alpha_n: {}",
                beta_n, alpha_n);
        }
    }

    double adaptive_dhat_ratio() const { return m_adaptive_dhat_ratio; }

    void set_adaptive_dhat_ratio(const double adaptive_dhat_ratio)
    {
        m_adaptive_dhat_ratio = adaptive_dhat_ratio;
    }

    double dhat = 1;
    double alpha_t = 1;
    double beta_t = 0;
    double alpha_n = 0.1;
    double beta_n = 0;
    int r = 2;

private:
    double m_adaptive_dhat_ratio = 0.5;
};

} // namespace ipc