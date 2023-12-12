#include <ipc/utils/quadrature.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/barrier/barrier.hpp>
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
        VectorMax3d tangent = (e1 - e0).normalized() / 1.0;
        VectorMax3d sample = e0 + (e1 - e0) * uv;
        VectorMax6d grad = point_point_distance_gradient(p, sample); // grad of dist^2 actually

        const double dist = sqrt(point_point_distance(p, sample));
        const double deriv = tangent.dot(grad.head(tangent.size())) / (2 * dist);
        const double dist_barrier = barrier(dist, dhat);
        
        double val = 0;
        for (double sign : {-1, 1})
        {
            val += dist_barrier * smooth_heaviside(deriv * sign + alpha);
        }

        return val;
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