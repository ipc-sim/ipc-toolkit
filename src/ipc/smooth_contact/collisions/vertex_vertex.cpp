#include "vertex_vertex.hpp"
#include <ipc/smooth_contact/pairs/smooth_point_point.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/AutodiffTypes.hpp>

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
        auto points = slice_positions<double, 6, 2>(positions);

        Vector<double, 2> direc = points.row(1) - points.row(0);
        const double dist = direc.norm();
        if (direc.norm() >= get_dhat())
            return false;
        direc /= dist;

        return smooth_point2_term_type(points.row(0), -direc, points.row(2), points.row(3), params.alpha, params.beta) &&
                smooth_point2_term_type(points.row(1), direc, points.row(4), points.row(5), params.alpha, params.beta);

    }

    template<class scalar>
    scalar SmoothVertexVertexCollision::evaluate_quadrature(const Vector12d& positions, ParameterType params) const
    {
        auto points = slice_positions<scalar, 6, 2>(positions);
        params.dhat = get_dhat();
        return smooth_point_point_potential_2d<scalar>(
            points.row(0), points.row(1), points.row(2), points.row(3), points.row(4), points.row(5), params);
    }

    double SmoothVertexVertexCollision::compute_distance(const Vector<double, -1, 18>& positions) const
    {
        auto points = slice_positions<double, 6, 2>(positions);

        return (points.row(0) - points.row(1)).squaredNorm();
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
        using Diff=ADGrad<12>;
        return evaluate_quadrature<Diff>(positions, params).getGradient();
    }

    MatrixMax18d SmoothVertexVertexCollision::hessian(
        const VectorMax18d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        DiffScalarBase::setVariableCount(12);
        using Diff=ADHessian<12>;
        return evaluate_quadrature<Diff>(positions, params).getHessian();
    }
}