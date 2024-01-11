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
    scalar SmoothEdgeEdgeCollision<dim_>::evaluate_quadrature(const VectorMax12d& positions, ParameterType params) const
    {
        std::array<VectorMax3<scalar>, 4> points = slice_positions<scalar>(positions, dim_);

        scalar val = scalar(0.);
        if constexpr (dim_ == 2)
        {
            if (params.n_quadrature > 1)
            {
                Eigen::VectorXd uv = Eigen::VectorXd::LinSpaced(params.n_quadrature, 0, 1);
                Vector2<scalar> p;
                std::array<Vector2<scalar>, 2> edges = {{points[1] - points[0], points[3] - points[2]}};
                const std::array<scalar, 2> lens = {{edges[0].norm(), edges[1].norm()}};

                for (const int e : {0, 1})
                {
                    const int ee = 1 - e;

                    if (dhat0 > 0 && dhat1 > 0)
                        params.eps = ee == 0 ? dhat0*dhat0 : dhat1*dhat1;
                    
                    scalar tmp = scalar(0.);
                    for (int i = 1; i < uv.size(); ++i)
                    {
                        for (int j = 1; j < uv.size(); ++j)
                        {
                            p = points[2*e] + scalar((uv(i) + uv(i-1)) / 2) * edges[e];
                            tmp += scalar(uv(i) - uv(i-1)) * smooth_point_edge_potential_single_point<scalar>(p, points[ee*2] + scalar(uv(j-1)) * edges[ee], points[ee*2] + scalar(uv(j)) * edges[ee], params);
                        }
                    }
                    val += tmp * lens[e];
                }
            }
            else
            {
                for (const int e : {0, 1})
                {
                    const int ee = 1 - e;

                    if (dhat0 > 0 && dhat1 > 0)
                        params.eps = ee == 0 ? dhat0*dhat0 : dhat1*dhat1;

                    const scalar len = (points[e * 2 + 1] - points[e * 2 + 0]).norm();
                    for (const int i : {0, 1})
                    {
                        const int p_id = e * 2 + i;
                        if (vertices[p_id] == vertices[ee * 2 + 0] ||
                            vertices[p_id] == vertices[ee * 2 + 1])
                            continue;
                        
                        val += (len / scalar(2.)) * smooth_point_edge_potential_single_point<scalar>(points[p_id], points[ee * 2 + 0], points[ee * 2 + 1], params);
                    }
                }
            }
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
    // template class SmoothEdgeEdgeCollision<3>;
}