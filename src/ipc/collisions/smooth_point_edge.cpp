#include <ipc/utils/quadrature.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/math.hpp>
#include <ipc/utils/AutodiffTypes.hpp>

const static int edge_quadrature_order = 1;

namespace ipc {
    template <typename scalar>
    scalar smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const scalar &uv,
        const double &dhat,
        const double &alpha)
    {
        VectorMax3<scalar> tangent = (e1 - e0).normalized();
        VectorMax3<scalar> sample = e0 + (e1 - e0) * uv;

        const scalar dist = (p - sample).norm();
        const scalar Phi_sqrt = tangent.dot((p - sample) / dist);

        return inv_barrier(dist, dhat, 1.) * cubic_spline(Phi_sqrt * 2. / alpha);
    }

    double smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3<double>>& p,
        const Eigen::Ref<const VectorMax3<double>>& e0,
        const Eigen::Ref<const VectorMax3<double>>& e1,
        const double &uv,
        const double &dhat,
        const double &alpha);
    AutodiffScalarGrad smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad>>& e1,
        const AutodiffScalarGrad &uv,
        const double &dhat,
        const double &alpha);
    AutodiffScalarHessian smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian>>& e1,
        const AutodiffScalarHessian &uv,
        const double &dhat,
        const double &alpha);

    template <typename scalar>
    scalar smooth_point_edge_potential(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const double &dhat,
        const double &alpha,
        const int &N)
    {
        Eigen::VectorXd pts, weights;
        ipc::line_quadrature(edge_quadrature_order, N, pts, weights);

        scalar val = 0;
        for (int i = 0; i < pts.size(); i++)
        {
            val += smooth_point_edge_potential_pointwise(p, e0, e1, pts(i), dhat, alpha) * weights(i);
        }

        return val;
    }

    double smooth_point_edge_potential(
        const Eigen::Ref<const VectorMax3<double>>& p,
        const Eigen::Ref<const VectorMax3<double>>& e0,
        const Eigen::Ref<const VectorMax3<double>>& e1,
        const double &dhat,
        const double &alpha,
        const int &N);
    AutodiffScalarGrad smooth_point_edge_potential(
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad>>& e1,
        const double &dhat,
        const double &alpha,
        const int &N);
    AutodiffScalarHessian smooth_point_edge_potential(
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian>>& e1,
        const double &dhat,
        const double &alpha,
        const int &N);
}