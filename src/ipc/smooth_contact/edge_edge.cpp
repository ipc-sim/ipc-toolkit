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
    const std::array<double, 2> &dhats_,
    const Eigen::MatrixXd &V): SmoothCollision<6>(primitive0_, primitive1_, dhats_, mesh)
    {
        vertices = vertex_ids(mesh.edges(), mesh.faces());
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

                if (Phi < params.alpha && dist_sqr < pow(get_dhat(ee), 2))
                {
                    normal_types[p_id] = HEAVISIDE_TYPE::VARIANT;
                    skip = false;
                }
            }
        }

        return !skip;
    }

    template<class scalar>
    scalar SmoothEdgeEdgeCollision::evaluate_quadrature(const Vector8d& positions, ParameterType params) const
    {
        std::array<Vector2<scalar>, 4> points = slice_positions<scalar, 4, 2>(positions);

        params.eps = get_eps();
        scalar val = scalar(0.);
        {
            for (const int e : {0, 1})
            {
                const int ee = 1 - e;

                params.eps = pow(get_dhat(ee), 2);

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

    double SmoothEdgeEdgeCollision::compute_distance(const Vector<double, -1, 18>& positions) const
    {
        std::array<Vector2<double>, 4> points = slice_positions<double, 4, 2>(positions);
        double min_dist = std::numeric_limits<double>::max();
        if (normal_types[0] != HEAVISIDE_TYPE::ZERO)
            min_dist = std::min(min_dist, point_edge_distance(points[0], points[2], points[3]));
        if (normal_types[1] != HEAVISIDE_TYPE::ZERO)
            min_dist = std::min(min_dist, point_edge_distance(points[1], points[2], points[3]));
        if (normal_types[2] != HEAVISIDE_TYPE::ZERO)
            min_dist = std::min(min_dist, point_edge_distance(points[2], points[0], points[1]));
        if (normal_types[3] != HEAVISIDE_TYPE::ZERO)
            min_dist = std::min(min_dist, point_edge_distance(points[3], points[0], points[1]));
        
        return min_dist;
    }

    double SmoothEdgeEdgeCollision::operator()(
        const VectorMax18d& positions, 
        const ParameterType &params) const
    {
        return evaluate_quadrature<double>(positions, params);
    }

    VectorMax18d SmoothEdgeEdgeCollision::gradient(
        const VectorMax18d& positions, 
        const ParameterType &params) const
    {
        DiffScalarBase::setVariableCount(8);
        using Diff=AutodiffScalarGrad<8>;
        return evaluate_quadrature<Diff>(positions, params).getGradient();
    }

    MatrixMax18d SmoothEdgeEdgeCollision::hessian(
        const VectorMax18d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        DiffScalarBase::setVariableCount(8);
        using Diff=AutodiffScalarHessian<8>;
        return evaluate_quadrature<Diff>(positions, params).getHessian();
    }
}