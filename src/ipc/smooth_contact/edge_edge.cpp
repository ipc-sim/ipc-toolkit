#include "edge_edge.hpp"
#include "smooth_point_edge.hpp"
#include <ipc/utils/math.hpp>
#include <ipc/utils/quadrature.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/utils/AutodiffTypes.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include <algorithm>

namespace ipc {

    SmoothEdgeEdgeCollision::SmoothEdgeEdgeCollision(
    long primitive0_,
    long primitive1_,
    const CollisionMesh &mesh,
    const ParameterType &param,
    const Eigen::MatrixXd &V): SmoothEdgeEdgeCollision(primitive0_, primitive1_, mesh)
    {
        Vector8d positions = dof(V, mesh.edges(), mesh.faces());
        is_active_ = compute_types(positions, param);
    }

    bool SmoothEdgeEdgeCollision::compute_types(const Vector8d& positions, const ParameterType &params)
    {
        std::array<Vector<double, 2>, 4> points = slice_positions<double, 4, 2>(positions);
        bool skip = true;
        normal_types.fill(HEAVISIDE_TYPE::ZERO);
        for (const int e : {0, 1})
        {
            const int ee = 1 - e;
            const auto &e0 = points[ee * 2 + 0];
            const auto &e1 = points[ee * 2 + 1];
            
            Vector<double, 2> tangent = e1 - e0;
            const double len = tangent.norm();
            tangent = tangent / len;

            for (const int i : {0, 1})
            {
                const int p_id = e * 2 + i;
                const auto &p = points[p_id];

                if (vertices[p_id] == vertices[ee * 2 + 0] || vertices[p_id] == vertices[ee * 2 + 1])
                    continue;

                const Vector<double, 2> pos = p - e0;
                const double s = pos.dot(tangent) / len;
                const double L = L_ns(s);
                const double dist_sqr = (pos - (L * len) * tangent).squaredNorm();
                const double Phi = 1 - cross2<double>(pos, tangent) / sqrt(dist_sqr);

                if (Phi < params.alpha && dist_sqr < params.eps)
                {
                    normal_types[p_id] = HEAVISIDE_TYPE::VARIANT;
                    skip = false;
                }
            }
        }

        return (params.n_quadrature > 1) || !skip;
    }

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

                    if (normal_types[p_id] != HEAVISIDE_TYPE::ZERO)
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

    double SmoothEdgeEdgeCollision::compute_distance(const Vector<double, -1, 12>& positions) const
    {
        std::array<Vector2<double>, 4> points = slice_positions<double, 4, 2>(positions);
        double min_dist = std::numeric_limits<double>::max();
        if (points[0] != points[2] && points[0] != points[3]) 
            min_dist = std::min(min_dist, point_edge_distance(points[0], points[2], points[3]));
        if (points[1] != points[2] && points[1] != points[3]) 
            min_dist = std::min(min_dist, point_edge_distance(points[1], points[2], points[3]));
        if (points[2] != points[0] && points[2] != points[1]) 
            min_dist = std::min(min_dist, point_edge_distance(points[2], points[0], points[1]));
        if (points[3] != points[0] && points[3] != points[1]) 
            min_dist = std::min(min_dist, point_edge_distance(points[3], points[0], points[1]));
        
        return min_dist;
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