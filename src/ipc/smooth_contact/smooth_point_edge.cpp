#include <ipc/utils/quadrature.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/math.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
#include "smooth_point_edge.hpp"
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <stdexcept>
#include <iostream>

namespace ipc {

    template <typename scalar>
    scalar smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const double &uv,
        const ParameterType &params)
    {
        VectorMax3<scalar> tangent = (e1 - e0).normalized();
        VectorMax3<scalar> sample = e0 + (e1 - e0) * scalar(uv);

        const scalar dist_sqr = (p - sample).squaredNorm();
        const scalar Phi = intpow(tangent.dot((p - sample)), 2) / dist_sqr;

        return inv_barrier(dist_sqr, params.eps, params.r) * cubic_spline(Phi * (2. / params.alpha));
    }

    template double smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3<double>>& p,
        const Eigen::Ref<const VectorMax3<double>>& e0,
        const Eigen::Ref<const VectorMax3<double>>& e1,
        const double &uv,
        const ParameterType &params);
    template AutodiffScalarGrad<12> smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& e1,
        const double &uv,
        const ParameterType &params);
    template AutodiffScalarHessian<12> smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& e1,
        const double &uv,
        const ParameterType &params);

    template <typename scalar>
    scalar smooth_point_edge_potential_quadrature(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const ParameterType &params)
    {
        Eigen::VectorXd pts, weights;
        ipc::line_quadrature(params.n_quadrature, pts, weights);

        scalar val(0);
        for (int i = 0; i < pts.size(); i++)
        {
            val += smooth_point_edge_potential_pointwise<scalar>(p, e0, e1, pts(i), params) * weights(i);
        }

        return val;
    }

    template double smooth_point_edge_potential_quadrature(
        const Eigen::Ref<const VectorMax3<double>>& p,
        const Eigen::Ref<const VectorMax3<double>>& e0,
        const Eigen::Ref<const VectorMax3<double>>& e1,
        const ParameterType &params);
    template AutodiffScalarGrad<12> smooth_point_edge_potential_quadrature(
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& e1,
        const ParameterType &params);
    template AutodiffScalarHessian<12> smooth_point_edge_potential_quadrature(
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& e1,
        const ParameterType &params);

    template <typename scalar>
    Vector<scalar, -1, -1> smooth_point_edge_potentials_quadrature(
        const Eigen::Ref<Eigen::Matrix<scalar, -1, -1, Eigen::RowMajor, -1, 3>>& points,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const ParameterType &params)
    {
        Vector<scalar, -1, -1> out;
        out.setZero(points.rows());

        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), points.rows()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i < r.end(); i++) {
                    out(i) = smooth_point_edge_potential_quadrature<scalar>(points.row(i).transpose(), e0, e1, params);
                }
            });
        
        return out;
    }

    template Vector<double, -1, -1> smooth_point_edge_potentials_quadrature(
        const Eigen::Ref<Eigen::Matrix<double, -1, -1, Eigen::RowMajor, -1, 3>>& points,
        const Eigen::Ref<const VectorMax3<double>>& e0,
        const Eigen::Ref<const VectorMax3<double>>& e1,
        const ParameterType &params);

    template <typename scalar>
    scalar smooth_point_edge_potential_single_point(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const ParameterType &params)
    {
        VectorMax3<scalar> tangent = e1 - e0;
        const scalar len = tangent.norm();
        tangent = tangent / len;

        const VectorMax3<scalar> pos = p - e0;
        const scalar s = pos.dot(tangent) / len;
        const scalar L = (params.a > 0) ? L_s(s, params.a) : L_ns(s);
        const VectorMax3<scalar> sample = e0 + (s - L) * len * tangent;
        const scalar Phi = intpow((sample - p).normalized().dot(tangent), 2);
        const scalar dist_sqr = intpow(cross2<scalar>(pos, tangent), 2) + intpow(len * L, 2);

        if (Phi > params.alpha)
            return scalar(0.);
        if (dist_sqr > params.eps)
            return scalar(0.);

        return len * cubic_spline(Phi * (2. / params.alpha)) * inv_barrier(dist_sqr, params.eps, params.r);
    }

    template double smooth_point_edge_potential_single_point(
        const Eigen::Ref<const VectorMax3<double>>& p,
        const Eigen::Ref<const VectorMax3<double>>& e0,
        const Eigen::Ref<const VectorMax3<double>>& e1,
        const ParameterType &params);
    template AutodiffScalarGrad<12> smooth_point_edge_potential_single_point(
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<12>>>& e1,
        const ParameterType &params);
    template AutodiffScalarHessian<12> smooth_point_edge_potential_single_point(
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<12>>>& e1,
        const ParameterType &params);
}