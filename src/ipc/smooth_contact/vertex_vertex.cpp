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
                else
                    vertices[2*j+3] = mesh.edges()(i, 0);
            }

        Vector12d positions = dof(V, mesh.edges(), mesh.faces());
        is_active_ = compute_types(positions, param);
    }

    bool SmoothVertexVertexCollision::compute_types(const Vector12d& positions, const ParameterType &params)
    {
        std::array<Vector<double, 2>, 6> points = slice_positions<double, 6, 2>(positions);

        if ((points[1] - points[0]).norm() >= std::max(get_dhat(0), get_dhat(1)))
            return false;

        bool all_skip = true;
        for (int j : {0, 1})
        {
            bool skip = false;
            for (int v : {0, 1})
            {
                const double val = (v == 1 ? -1 : 1) * (points[j] - points[1-j]).normalized().dot((points[2 + j * 2 + v] - points[j]).normalized());
                if (val < -params.alpha)
                    skip = true;
            }
            if (!skip)
                all_skip = false;
        }

        return !all_skip;
    }

    template<class scalar>
    scalar SmoothVertexVertexCollision::evaluate_quadrature(const Vector12d& positions, ParameterType params) const
    {
        std::array<Vector2<scalar>, 6> points = slice_positions<scalar, 6, 2>(positions);
        return smooth_point_point_potential<scalar>(
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