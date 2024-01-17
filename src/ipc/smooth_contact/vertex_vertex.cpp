#include "vertex_vertex.hpp"
#include "smooth_point_point.hpp"
#include <ipc/utils/math.hpp>
#include <ipc/utils/quadrature.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/AutodiffTypes.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include <algorithm>

namespace ipc {

    SmoothVertexVertexCollision::SmoothVertexVertexCollision(
    long primitive0_,
    long primitive1_,
    const CollisionMesh &mesh,
    const ParameterType &param,
    const std::array<double, 2> &dhats_,
    const Eigen::MatrixXd &V): SmoothCollision<6>(primitive0_, primitive1_, dhats_, mesh)
    {
        vertices[0] = primitive0_;
        vertices[1] = primitive1_;
        for (int j : {0, 1})
            for (long i : mesh.vertex_edge_adjacencies()[vertices[j]])
            {
                if (mesh.edges()(i, 0) == vertices[j])
                    vertices[2*j+2] = mesh.edges()(i, 1);
                else if (mesh.edges()(i, 1) == vertices[j])
                    vertices[2*j+3] = mesh.edges()(i, 0);
                else
                    throw std::runtime_error("Invalid edge-vertex adjacency!");
            }

        Vector12d positions = dof(V, mesh.edges(), mesh.faces());
        is_active_ = compute_types(positions, param);
    }

    bool SmoothVertexVertexCollision::compute_types(const Vector12d& positions, const ParameterType &params)
    {
        std::array<Vector<double, 2>, 6> points = slice_positions<double, 6, 2>(positions);

        Vector<double, 2> direc = points[1] - points[0];
        Vector<double, 2> ta0 = points[2] - points[0], ta1 = points[0] - points[3];
        Vector<double, 2> tb0 = points[4] - points[1], tb1 = points[1] - points[5];
        ta0.normalize(); ta1.normalize();
        tb0.normalize(); tb1.normalize();

        if (direc.norm() >= std::max(get_dhat(0), get_dhat(1)))
            return false;
        direc.normalize();
        
        // tangent term
        {
            bool A = -direc.dot(ta0) <= -params.alpha || -direc.dot(-ta1) <= -params.alpha;
            bool B = direc.dot(tb0) <= -params.alpha || direc.dot(-tb1) <= -params.alpha;
            if (A && B)
                return false;
        }

        // normal term
        {
            bool A = cross2<double>(direc, ta0) <= -params.alpha && cross2<double>(direc, ta1) <= -params.alpha;
            bool B = cross2<double>(direc, tb0) >= params.alpha && cross2<double>(direc, tb1) >= params.alpha;
            if (A || B)
                return false;
        }

        return true;
    }

    template<class scalar>
    scalar SmoothVertexVertexCollision::evaluate_quadrature(const Vector12d& positions, ParameterType params) const
    {
        std::array<Vector2<scalar>, 6> points = slice_positions<scalar, 6, 2>(positions);
        return smooth_point_point_potential_2d<scalar>(
            points[0], points[1], points[2], points[3], points[4], points[5], params, dhats);
    }

    double SmoothVertexVertexCollision::compute_distance(const Vector<double, -1, 18>& positions) const
    {
        std::array<Vector2<double>, 6> points = slice_positions<double, 6, 2>(positions);

        return (points[0] - points[1]).squaredNorm();
    }

    double SmoothVertexVertexCollision::operator()(
        const VectorMax18d& positions, 
        const ParameterType &params) const
    {
        return evaluate_quadrature<double>(positions, params);
    }

    VectorMax18d SmoothVertexVertexCollision::gradient(
        const VectorMax18d& positions, 
        const ParameterType &params) const
    {
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarGrad<12>;
        return evaluate_quadrature<Diff>(positions, params).getGradient();
    }

    MatrixMax18d SmoothVertexVertexCollision::hessian(
        const VectorMax18d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarHessian<12>;
        return evaluate_quadrature<Diff>(positions, params).getHessian();
    }
}