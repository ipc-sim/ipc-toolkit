#include "edge_vertex.hpp"
#include <ipc/smooth_contact/pairs/smooth_point_edge.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/utils/quadrature.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
#include <iostream>

DECLARE_DIFFSCALAR_BASE();

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include <algorithm>

namespace ipc {

SmoothEdgeVertexCollision::SmoothEdgeVertexCollision(
    long primitive0_,
    long primitive1_,
    const CollisionMesh& mesh,
    const ParameterType& param,
    const double& dhat,
    const Eigen::MatrixXd& V)
    : SmoothCollision<6>(primitive0_, primitive1_, dhat, mesh)
{
    vertices[0] = primitive1;
    vertices[1] = mesh.edges()(primitive0, 0);
    vertices[2] = mesh.edges()(primitive0, 1);
    for (long i : mesh.vertex_edge_adjacencies()[vertices[0]]) {
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

bool SmoothEdgeVertexCollision::compute_types(
    const Vector10d& positions, const ParameterType& params)
{
    auto points = slice_positions<double, 5, 2>(positions);

    auto dtype =
        point_edge_distance_type(points.row(0), points.row(1), points.row(2));
    if (dtype != PointEdgeDistanceType::P_E)
        return false;

    Vector<double, 2> direc =
        PointEdgeDistance<double, 2>::point_edge_closest_point_direction(
            points.row(0), points.row(1), points.row(2), dtype);
    const double dist = direc.norm();
    if (dist >= get_dhat())
        return false;
    direc /= dist;

    // edge term
    if (Math<double>::cross2(
            direc, (points.row(2) - points.row(1)).normalized())
        < 0) {
        // std::cout << "edge term zero\n";
        if (dist < 1e-10)
            std::cout << direc.transpose() << ", " << points.row(1).transpose()
                      << ", " << points.row(2).transpose() << "\n";
        return false;
    }

    // point term
    if (!smooth_point2_term_type(
            points.row(0), -direc, points.row(3), points.row(4), params.alpha,
            params.beta)) {
        // std::cout << "point term zero\n";
        if (dist < 1e-10)
            std::cout << direc.transpose() << ", " << points.row(1).transpose()
                      << ", " << points.row(2).transpose() << "\n";
        return false;
    }

    return true;
}

template <class scalar>
scalar SmoothEdgeVertexCollision::evaluate_quadrature(
    const Vector10d& positions, ParameterType params) const
{
    auto points = slice_positions<scalar, 5, 2>(positions);
    params.dhat = get_dhat();
    return smooth_point_edge_potential_single_point<scalar>(
        points.row(0), points.row(1), points.row(2), points.row(3),
        points.row(4), params);
}

double SmoothEdgeVertexCollision::compute_distance(
    const Vector<double, -1, 18>& positions) const
{
    auto points = slice_positions<double, 5, 2>(positions);

    return point_edge_distance(
        points.row(0), points.row(1), points.row(2),
        PointEdgeDistanceType::AUTO);
}

double SmoothEdgeVertexCollision::operator()(
    const VectorMax18d& positions, const ParameterType& params) const
{
    return evaluate_quadrature<double>(positions, params);
}

VectorMax18d SmoothEdgeVertexCollision::gradient(
    const VectorMax18d& positions, const ParameterType& params) const
{
    DiffScalarBase::setVariableCount(10);
    using Diff = ADGrad<10>;
    return evaluate_quadrature<Diff>(positions, params).getGradient();
}

MatrixMax18d SmoothEdgeVertexCollision::hessian(
    const VectorMax18d& positions, const ParameterType& params) const
{
    DiffScalarBase::setVariableCount(10);
    using Diff = ADHessian<10>;
    return evaluate_quadrature<Diff>(positions, params).getHessian();
}
} // namespace ipc