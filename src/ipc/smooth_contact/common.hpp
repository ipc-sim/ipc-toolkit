#pragma once
#include <cmath>

// #define DERIVATIVES_WITH_AUTODIFF

namespace ipc {

constexpr static int n_vert_neighbors_2d = 3;
constexpr static int n_edge_neighbors_2d = 2;
constexpr static int max_vert_2d =
    2 * std::max(n_vert_neighbors_2d, n_edge_neighbors_2d);
constexpr static int n_vert_neighbors_3d = 15;
constexpr static int n_edge_neighbors_3d = 4;
constexpr static int n_face_neighbors_3d = 3;
constexpr static int max_vert_3d = 24;

struct ParameterType {
    ParameterType(
        const double& _dhat,
        const double& _alpha,
        const double& _r,
        const double& _beta)
        : dhat(_dhat)
        , alpha(_alpha)
        , r(_r)
        , beta(_beta)
    {
        if (!(r > 0) || !(dhat > 0) || !(alpha > 0) || !(beta + alpha > 1e-8))
            logger().error(
                "Wrong parameters for smooth contact! dhat {} alpha {} r {} beta {}",
                dhat, alpha, r, beta);
    }
    ParameterType() = delete;

    void set_adaptive_dhat_ratio(const double adaptive_dhat_ratio_)
    {
        adaptive_dhat_ratio = adaptive_dhat_ratio_;
    }
    double get_adaptive_dhat_ratio() const { return adaptive_dhat_ratio; }

    double dhat;
    const double alpha;
    const double r;
    const double beta;

private:
    double adaptive_dhat_ratio = 0.5;
};

} // namespace ipc