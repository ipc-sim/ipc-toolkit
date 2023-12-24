#include <ipc/utils/quadrature.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/math.hpp>

const static int edge_quadrature_order = 1;

namespace ipc {
    double smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3d>& p,
        const Eigen::Ref<const VectorMax3d>& e0,
        const Eigen::Ref<const VectorMax3d>& e1,
        const double &uv,
        const double &dhat,
        const double &alpha)
    {
        VectorMax3d tangent = (e1 - e0).normalized();
        VectorMax3d sample = e0 + (e1 - e0) * uv;

        const double dist = (p - sample).norm();
        const double Phi_sqrt = tangent.dot((p - sample).normalized());

        return inv_barrier(dist, dhat, 2) * cubic_spline(Phi_sqrt * 2. / alpha);
    }

    double smooth_point_edge_potential(
        const Eigen::Ref<const VectorMax3d>& p,
        const Eigen::Ref<const VectorMax3d>& e0,
        const Eigen::Ref<const VectorMax3d>& e1,
        const double &dhat,
        const double &alpha,
        const int &N)
    {
        Eigen::VectorXd pts, weights;
        ipc::line_quadrature(edge_quadrature_order, N, pts, weights);

        double val = 0;
        for (int i = 0; i < pts.size(); i++)
        {
            val += smooth_point_edge_potential_pointwise(p, e0, e1, pts(i), dhat, alpha) * weights(i);
        }

        return val;
    }
}