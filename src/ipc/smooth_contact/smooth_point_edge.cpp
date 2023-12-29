#include <ipc/utils/quadrature.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/math.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
#include "smooth_point_edge.hpp"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace ipc {
    double SmoothEdgeVertexConstraint::compute_potential(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const SmoothParameter param) const
    {
        return smooth_point_edge_potential_single_point<double>(vertices.row(vertex_id), vertices.row(edges(edge_id, 0)), vertices.row(edges(edge_id, 1)), param.eps, param.alpha, param.r);
    }

    VectorMax12d SmoothEdgeVertexConstraint::compute_potential_gradient(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const SmoothParameter param) const
    {
        const int dim = vertices.cols();
        DiffScalarBase::setVariableCount(9);
        using Diff=AutodiffScalarGrad<9>;
        VectorMax3<Diff> p, e0, e1;
        p.setZero(dim);
        e0.setZero(dim);
        e1.setZero(dim);
        for (int d = 0; d < dim; d++)
        {
            p(d) = Diff(d, vertices(vertex_id, d));
            e0(d) = Diff(3 + d, vertices(edges(edge_id, 0), d));
            e1(d) = Diff(6 + d, vertices(edges(edge_id, 1), d));
        }

        const auto val = smooth_point_edge_potential_single_point<Diff>(p, e0, e1, param.eps, param.alpha, param.r);

        VectorMax12d grad;
        grad.setZero(3 * dim);
        for (int d = 0; d < dim; d++)
            for (int i = 0; i < 3; i++)
                grad(d + dim * i) = val.getGradient()(d + 3 * i);
        
        return grad;
    }

    MatrixMax12d SmoothEdgeVertexConstraint::compute_potential_hessian(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const SmoothParameter param,
        const bool project_hessian_to_psd) const
    {
        const int dim = vertices.cols();
        DiffScalarBase::setVariableCount(9);
        using Diff=AutodiffScalarHessian<9>;
        VectorMax3<Diff> p, e0, e1;
        p.setZero(dim);
        e0.setZero(dim);
        e1.setZero(dim);
        for (int d = 0; d < dim; d++)
        {
            p(d) = Diff(d, vertices(vertex_id, d));
            e0(d) = Diff(3 + d, vertices(edges(edge_id, 0), d));
            e1(d) = Diff(6 + d, vertices(edges(edge_id, 1), d));
        }

        const auto val = smooth_point_edge_potential_single_point<Diff>(p, e0, e1, param.eps, param.alpha, param.r);

        MatrixMax12d hess;
        hess.setZero(3 * dim, 3 * dim);
        for (int d1 = 0; d1 < dim; d1++)
            for (int i1 = 0; i1 < 3; i1++)
                for (int d2 = 0; d2 < dim; d2++)
                    for (int i2 = 0; i2 < 3; i2++)
                        hess(d1 + dim * i1, d2 + dim * i2) = val.getHessian()(d1 + 3 * i1, d2 + 3 * i2);
        
        return hess;
    }

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
    template AutodiffScalarGrad<9> smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<9>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<9>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<9>>>& e1,
        const double &uv,
        const double &dhat,
        const double &alpha);
    template AutodiffScalarHessian<9> smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<9>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<9>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<9>>>& e1,
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
    template AutodiffScalarGrad<9> smooth_point_edge_potential_quadrature(
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<9>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<9>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<9>>>& e1,
        const double &dhat,
        const double &alpha,
        const int &N);
    template AutodiffScalarHessian<9> smooth_point_edge_potential_quadrature(
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<9>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<9>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<9>>>& e1,
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
        const VectorMax3<scalar> tangent = e1 - e0;
        const scalar len = tangent.norm();

        assert(p.size() == 2); // normal only considered for 2d
        VectorMax3<scalar> normal = Vector2<scalar>::Zero();
        normal(0) = -tangent(1);
        normal(1) = tangent(0);

        const scalar s = (p - e0).dot(tangent) / pow(len, 2);
        const scalar L = L_s(s, 2);
        const VectorMax3<scalar> sample = e0 + (s - L) * tangent;
        const scalar Phi = pow((sample - p).dot(tangent) / len, 2);
        const scalar dist_sqr = pow((p - e0).dot(normal) / len, 2) + pow(len * L, 2);

        return len * cubic_spline(Phi * (2 / alpha)) * inv_barrier(dist_sqr, dhat, r);
    }

    template double smooth_point_edge_potential_single_point(
        const Eigen::Ref<const VectorMax3<double>>& p,
        const Eigen::Ref<const VectorMax3<double>>& e0,
        const Eigen::Ref<const VectorMax3<double>>& e1,
        const double &dhat,
        const double &alpha,
        const double &r);
    template AutodiffScalarGrad<9> smooth_point_edge_potential_single_point(
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<9>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<9>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarGrad<9>>>& e1,
        const double &dhat,
        const double &alpha,
        const double &r);
    template AutodiffScalarHessian<9> smooth_point_edge_potential_single_point(
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<9>>>& p,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<9>>>& e0,
        const Eigen::Ref<const VectorMax3<AutodiffScalarHessian<9>>>& e1,
        const double &dhat,
        const double &alpha,
        const double &r);

}