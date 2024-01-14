#include "edge_edge.hpp"
#include "smooth_point_edge.hpp"
#include <ipc/utils/math.hpp>
#include <ipc/utils/quadrature.hpp>
#include <ipc/utils/AutodiffTypes.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace ipc {

    template<class scalar>
    scalar SmoothEdgeEdgeCollision::evaluate_quadrature(const Vector8d& positions, ParameterType params) const
    {
        std::array<Vector2<scalar>, 4> points = slice_positions<scalar, 4, 2>(positions);

        scalar val = scalar(0.);
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

        return val;
    }

    Vector12d SmoothEdgeEdgeCollision::positions_to_3d(const Vector8d& positions) const
    {
        Vector12d positions_full;
        positions_full.setZero();
        for (int i = 0; i < 4; i++)
            positions_full.segment<dim>(i * 3) = positions.segment<dim>(i * dim);
        return positions_full;
    }

    double SmoothEdgeEdgeCollision::operator()(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        return evaluate_quadrature<double>(positions, params);
    }

    VectorMax12d SmoothEdgeEdgeCollision::gradient(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarGrad<12>;
        return evaluate_quadrature<Diff>(positions, params).getGradient().head(4*dim);
    }

    MatrixMax12d SmoothEdgeEdgeCollision::hessian(
        const VectorMax12d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarHessian<12>;
        return evaluate_quadrature<Diff>(positions, params).getHessian().topLeftCorner(4*dim, 4*dim);
    }
}