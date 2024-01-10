#include "edge_edge.hpp"
#include "smooth_point_edge.hpp"
#include <ipc/utils/quadrature.hpp>
#include <ipc/utils/AutodiffTypes.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace ipc {
    namespace {
        template <class T>
        std::array<VectorMax3<T>, 4> slice_positions(const VectorMax12d &positions, const int &dim_)
        {
            std::array<VectorMax3<T>, 4> points;
            points.fill(VectorMax3<T>::Zero(dim_));
            
            for (int i = 0, id = 0; i < 4; i++)
                for (int d = 0; d < dim_; d++, id++)
                    if constexpr (std::is_same<T, double>::value)
                        points[i](d) = positions(id);
                    else
                        points[i](d) = T(id, positions(id));

            return points;
        }
    }

    template<int dim_> template<class scalar>
    scalar SmoothEdgeEdgeCollision<dim_>::evaluate_quadrature(const VectorMax12d& positions, const ParameterType &params) const
    {
        auto [e00, e01, e10, e11] = slice_positions<scalar>(positions, dim_);

        scalar val = scalar(0.);
        if constexpr (dim_ == 2)
        {
            Eigen::VectorXd uv, w;
            line_quadrature(params.n_quadrature, uv, w);
            Vector2<scalar> p;
            const Vector2<scalar> edge0 = e01 - e00;
            const Vector2<scalar> edge1 = e11 - e10;
            const scalar len0 = edge0.norm(), len1 = edge1.norm();

            assert(uv(uv.size()-1) - uv(0) > 1 - 1e-10);
            
            scalar val0 = scalar(0.);
            scalar val1 = scalar(0.);
            for (int i = 0; i < uv.size(); ++i)
            {
                for (int j = 1; j < uv.size(); ++j)
                {
                    p = e00 + scalar(uv(i)) * edge0;
                    val0 += scalar(w(i)) * smooth_point_edge_potential_single_point<scalar>(p, e10 + scalar(uv(j-1)) * edge1, e10 + scalar(uv(j)) * edge1, params);

                    p = e10 + scalar(uv(i)) * edge1;
                    val1 += scalar(w(i)) * smooth_point_edge_potential_single_point<scalar>(p, e00 + scalar(uv(j-1)) * edge0, e00 + scalar(uv(j)) * edge0, params);
                }
            }
            val = val0 * len0 + val1 * len1;
        }

        return val;
    }

    template <int dim_>
    Vector12d SmoothEdgeEdgeCollision<dim_>::positions_to_3d(const VectorMax12d& positions) const
    {
        if constexpr (dim_ == 3)
            return positions;
        
        assert(positions.size() == 8);
        Vector12d positions_full;
        positions_full.setZero();
        for (int i = 0; i < 4; i++)
            positions_full.segment<dim_>(i * 3) = positions.segment<dim_>(i * dim_);
        return positions_full;
    }

    template <int dim_>
    double SmoothEdgeEdgeCollision<dim_>::compute_distance(const VectorMax12d& positions) const
    {
        if constexpr (dim_ == 3)
            return EdgeEdgeCandidate::compute_distance(positions);
        else
            return EdgeEdgeCandidate::compute_distance(positions_to_3d(positions));
    }

    template <int dim_>
    VectorMax12d SmoothEdgeEdgeCollision<dim_>::compute_distance_gradient(
        const VectorMax12d& positions) const
    {
        if constexpr (dim_ == 3)
            return EdgeEdgeCandidate::compute_distance_gradient(positions);
        else
            return EdgeEdgeCandidate::compute_distance_gradient(positions_to_3d(positions));
    }

    template <int dim_> MatrixMax12d
    SmoothEdgeEdgeCollision<dim_>::compute_distance_hessian(const VectorMax12d& positions) const
    {
        if constexpr (dim_ == 3)
            return EdgeEdgeCandidate::compute_distance_hessian(positions);
        else
            return EdgeEdgeCandidate::compute_distance_hessian(positions_to_3d(positions));
    }

    template <int dim_>
    double SmoothEdgeEdgeCollision<dim_>::operator()(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        return evaluate_quadrature<double>(positions, params);
    }

    template <int dim_>
    VectorMax12d SmoothEdgeEdgeCollision<dim_>::gradient(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarGrad<12>;
        return evaluate_quadrature<Diff>(positions, params).getGradient().head(4*dim_);
    }

    template <int dim_>
    MatrixMax12d SmoothEdgeEdgeCollision<dim_>::hessian(
        const VectorMax12d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarHessian<12>;
        return evaluate_quadrature<Diff>(positions, params).getHessian().topLeftCorner(4*dim_, 4*dim_);
    }

    template class SmoothEdgeEdgeCollision<2>;
    template class SmoothEdgeEdgeCollision<3>;
}