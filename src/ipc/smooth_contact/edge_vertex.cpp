#include "edge_vertex.hpp"
#include "smooth_point_edge.hpp"
#include <ipc/utils/math.hpp>
#include <ipc/utils/quadrature.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
#include <ipc/utils/distance_autodiff.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include <algorithm>

namespace ipc {

    SmoothEdgeVertexCollision::SmoothEdgeVertexCollision(
    long primitive0_,
    long primitive1_,
    const CollisionMesh &mesh,
    const ParameterType &param,
    const std::array<double, 2> &dhats_,
    const Eigen::MatrixXd &V): SmoothCollision<6>(primitive0_, primitive1_, dhats_, mesh)
    {
        vertices[0] = primitive1;
        vertices[1] = mesh.edges()(primitive0, 0);
        vertices[2] = mesh.edges()(primitive0, 1);
        for (long i : mesh.vertex_edge_adjacencies()[vertices[0]])
        {
            if (mesh.edges()(i, 0) == vertices[0])
                vertices[3] = mesh.edges()(i, 1);
            else if (mesh.edges()(i, 1) == vertices[0])
                vertices[4] = mesh.edges()(i, 0);
            else
                throw std::runtime_error("Invalid edge-vertex adjacency!");
        }

        Vector10d positions = dof(V, mesh.edges(), mesh.faces());
        is_active_ = compute_types(positions, param);
    }

    bool SmoothEdgeVertexCollision::compute_types(const Vector10d& positions, const ParameterType &params)
    {
        std::array<Vector<double, 2>, 5> points = slice_positions<double, 5, 2>(positions);

        Vector<double, 2> direc = point_edge_closest_point_direction<double>(points[0], points[1], points[2], PointEdgeDistanceType::AUTO);
        Vector<double, 2> e = (points[2] - points[1]).normalized();
        Vector<double, 2> t0 = (points[3] - points[0]).normalized(), t1 = (points[0] - points[4]).normalized();

        const double dist = direc.norm();
        if (dist >= std::min(get_dhat(0), get_dhat(1)))
            return false;
        direc /= dist;
        
        // tangent term
        {
            const double Phi = 1 - cross2<double>(points[0] - points[1], e) / dist;
            if (t0.dot(direc) <= -params.alpha || t1.dot(direc) <= -params.alpha || Phi >= params.alpha)
                return false;
        }

        // normal term
        {
            if (cross2<double>(direc, t0) <= -params.alpha && cross2<double>(direc, t1) <= -params.alpha)
                return false;
        }

        return true;
    }

    template<class scalar>
    scalar SmoothEdgeVertexCollision::evaluate_quadrature(const Vector10d& positions, ParameterType params) const
    {
        std::array<Vector2<scalar>, 5> points = slice_positions<scalar, 5, 2>(positions);
        params.eps = pow(std::min(get_dhat(0), get_dhat(1)), 2);
        return smooth_point_edge_potential_single_point<scalar>(
            points[0], points[1], points[2], points[3], points[4], params);
    }

    double SmoothEdgeVertexCollision::compute_distance(const Vector<double, -1, 18>& positions) const
    {
        std::array<Vector2<double>, 5> points = slice_positions<double, 5, 2>(positions);

        return point_edge_sqr_distance<double>(points[0], points[1], points[2], PointEdgeDistanceType::AUTO); 
    }

    double SmoothEdgeVertexCollision::operator()(
        const VectorMax18d& positions, 
        const ParameterType &params) const
    {
        return evaluate_quadrature<double>(positions, params);
    }

    VectorMax18d SmoothEdgeVertexCollision::gradient(
        const VectorMax18d& positions, 
        const ParameterType &params) const
    {
        DiffScalarBase::setVariableCount(10);
        using Diff=AutodiffScalarGrad<10>;
        return evaluate_quadrature<Diff>(positions, params).getGradient();
    }

    MatrixMax18d SmoothEdgeVertexCollision::hessian(
        const VectorMax18d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        DiffScalarBase::setVariableCount(10);
        using Diff=AutodiffScalarHessian<10>;
        return evaluate_quadrature<Diff>(positions, params).getHessian();
    }
}