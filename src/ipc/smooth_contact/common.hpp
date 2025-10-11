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
    static constexpr int VALUE = MAX_VERT_2D;
};
template <> class MaxVertices<3> {
public:
    static constexpr int VALUE = MAX_VERT_3D;
};

struct ParameterType {
    ParameterType(
        const double _dhat,
        const double _alpha_t,
        const double _beta_t,
        const int _r)
        : ParameterType(_dhat, _alpha_t, _beta_t, 0, 0.1, _r)
    {
    }
    ParameterType(
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
        if (!(r > 0) || !(dhat > 0) || !(abs(alpha_t) <= 1)
            || !(abs(alpha_n) <= 1) || !(abs(beta_t) <= 1)
            || !(abs(beta_n) <= 1) || !(beta_t + alpha_t > 1e-6)
            || !(beta_n + alpha_n > 1e-6)) {
            logger().error(
                "Wrong parameters for smooth contact! dhat {} alpha_t {} beta_t {} alpha_n {} beta_n {} r {}",
                dhat, alpha_t, beta_t, alpha_n, beta_n, r);
        }
    }
    ParameterType() { }

    double adaptive_dhat_ratio() const { return m_adaptive_dhat_ratio; }

    void set_adaptive_dhat_ratio(const double adaptive_dhat_ratio)
    {
        m_adaptive_dhat_ratio = adaptive_dhat_ratio;
    }

    double dhat = 1;
    double alpha_t = 1, beta_t = 0;
    double alpha_n = 0.1, beta_n = 0;
    int r = 2;

private:
    double m_adaptive_dhat_ratio = 0.5;
};

} // namespace ipc