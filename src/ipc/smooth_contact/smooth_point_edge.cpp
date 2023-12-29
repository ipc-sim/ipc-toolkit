#include <ipc/utils/quadrature.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/math.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
#include "smooth_point_edge.hpp"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace ipc {

    template <typename scalar>
    scalar smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const double &uv,
        const double &dhat,
        const double &alpha)
    {
        VectorMax3<scalar> tangent = (e1 - e0).normalized();
        VectorMax3<scalar> sample = e0 + (e1 - e0) * scalar(uv);

        const scalar dist = (p - sample).norm();
        const scalar Phi_sqrt = tangent.dot((p - sample) / dist);

        return inv_barrier(dist, dhat, 2.) * cubic_spline(Phi_sqrt * 2. / alpha);
    }

    template double smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3<double>>& p,
        const Eigen::Ref<const VectorMax3<double>>& e0,
        const Eigen::Ref<const VectorMax3<double>>& e1,
        const double &uv,
        const double &dhat,
        const double &alpha);
    template AutodiffScalarGrad<12> smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& e1,
        const double &uv,
        const double &dhat,
        const double &alpha);
    template AutodiffScalarHessian<12> smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& e1,
        const double &uv,
        const double &dhat,
        const double &alpha);

    template <typename scalar>
    scalar smooth_point_edge_potential_quadrature(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const double &dhat,
        const double &alpha,
        const int &N)
    {
        Eigen::VectorXd pts, weights;
        ipc::line_quadrature(N, pts, weights);

        scalar val(0);
        for (int i = 0; i < pts.size(); i++)
        {
            val += smooth_point_edge_potential_pointwise<scalar>(p, e0, e1, pts(i), dhat, alpha) * weights(i);
        }

        return val;
    }

    template double smooth_point_edge_potential_quadrature(
        const Eigen::Ref<const VectorMax3<double>>& p,
        const Eigen::Ref<const VectorMax3<double>>& e0,
        const Eigen::Ref<const VectorMax3<double>>& e1,
        const double &dhat,
        const double &alpha,
        const int &N);
    template AutodiffScalarGrad<12> smooth_point_edge_potential_quadrature(
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& e1,
        const double &dhat,
        const double &alpha,
        const int &N);
    template AutodiffScalarHessian<12> smooth_point_edge_potential_quadrature(
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& e1,
        const double &dhat,
        const double &alpha,
        const int &N);

    template <typename scalar>
    Vector<scalar, -1, -1> smooth_point_edge_potentials_quadrature(
        const Eigen::Ref<Eigen::Matrix<scalar, -1, -1, Eigen::RowMajor, -1, 3>>& points,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const double &dhat,
        const double &alpha,
        const int &N)
    {
        Vector<scalar, -1, -1> out;
        out.setZero(points.rows());

        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), points.rows()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i < r.end(); i++) {
                    out(i) = smooth_point_edge_potential_quadrature<scalar>(points.row(i).transpose(), e0, e1, dhat, alpha, N);
                }
            });
        
        return out;
    }

    template Vector<double, -1, -1> smooth_point_edge_potentials_quadrature(
        const Eigen::Ref<Eigen::Matrix<double, -1, -1, Eigen::RowMajor, -1, 3>>& points,
        const Eigen::Ref<const VectorMax3<double>>& e0,
        const Eigen::Ref<const VectorMax3<double>>& e1,
        const double &dhat,
        const double &alpha,
        const int &N);

    template <typename scalar>
    scalar smooth_point_edge_potential_single_point(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const double &dhat,
        const double &alpha,
        const double &r)
    {
        VectorMax3<scalar> tangent = e1 - e0;
        const scalar len = tangent.norm();
        tangent = tangent / len;

        const VectorMax3<scalar> pos = p - e0;
        const scalar s = pos.dot(tangent) / len;
        const scalar L = L_s(s, 2);
        const VectorMax3<scalar> sample = e0 + (s - L) * len * tangent;
        const scalar Phi = pow((sample - p).dot(tangent), 2);
        const scalar dist_sqr = pow(cross2<scalar>(pos, tangent), 2) + pow(len * L, 2);

        return len * cubic_spline(Phi * (2 / alpha)) * inv_barrier(dist_sqr, dhat, r);
    }

    template double smooth_point_edge_potential_single_point(
        const Eigen::Ref<const VectorMax3<double>>& p,
        const Eigen::Ref<const VectorMax3<double>>& e0,
        const Eigen::Ref<const VectorMax3<double>>& e1,
        const double &dhat,
        const double &alpha,
        const double &r);
    template AutodiffScalarGrad<12> smooth_point_edge_potential_single_point(
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& e1,
        const double &dhat,
        const double &alpha,
        const double &r);
    template AutodiffScalarHessian<12> smooth_point_edge_potential_single_point(
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& e1,
        const double &dhat,
        const double &alpha,
        const double &r);

}