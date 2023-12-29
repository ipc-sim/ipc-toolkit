#include "edge_edge.hpp"

#include <ipc/barrier/barrier.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>

namespace ipc {

EdgeEdgeCollision::EdgeEdgeCollision(
    const long _edge0_id,
    const long _edge1_id,
    const double _eps_x,
    const EdgeEdgeDistanceType _dtype)
    : EdgeEdgeCandidate(_edge0_id, _edge1_id)
    , eps_x(_eps_x)
    , dtype(_dtype)
{
}

EdgeEdgeCollision::EdgeEdgeCollision(
    const EdgeEdgeCandidate& candidate,
    const double _eps_x,
    const EdgeEdgeDistanceType _dtype)
    : EdgeEdgeCandidate(candidate)
    , eps_x(_eps_x)
    , dtype(_dtype)
{
}

EdgeEdgeCollision::EdgeEdgeCollision(
    const long _edge0_id,
    const long _edge1_id,
    const double _eps_x,
    const double _weight,
    const Eigen::SparseVector<double>& _weight_gradient,
    const EdgeEdgeDistanceType _dtype)
    : EdgeEdgeCandidate(_edge0_id, _edge1_id)
    , Collision(_weight, _weight_gradient)
    , eps_x(_eps_x)
    , dtype(_dtype)
{
}

double
EdgeEdgeCollision::mollifier_threshold(const VectorMax12d& rest_positions) const
{
    return edge_edge_mollifier_threshold(
        rest_positions.segment<3>(0), rest_positions.segment<3>(3),
        rest_positions.segment<3>(6), rest_positions.segment<3>(9));
}

double EdgeEdgeCollision::mollifier(const VectorMax12d& positions) const
{
    return mollifier(positions, eps_x);
}

double EdgeEdgeCollision::mollifier(
    const VectorMax12d& positions, const double _eps_x) const
{
    assert(positions.size() == 12);
    return edge_edge_mollifier(
        positions.segment<3>(0), positions.segment<3>(3),
        positions.segment<3>(6), positions.segment<3>(9), _eps_x);
}

VectorMax12d
EdgeEdgeCollision::mollifier_gradient(const VectorMax12d& positions) const
{
    return mollifier_gradient(positions, eps_x);
}

VectorMax12d EdgeEdgeCollision::mollifier_gradient(
    const VectorMax12d& positions, const double _eps_x) const
{
    assert(positions.size() == 12);
    return edge_edge_mollifier_gradient(
        positions.segment<3>(0), positions.segment<3>(3),
        positions.segment<3>(6), positions.segment<3>(9), _eps_x);
}

MatrixMax12d
EdgeEdgeCollision::mollifier_hessian(const VectorMax12d& positions) const
{
    return mollifier_hessian(positions, eps_x);
}

MatrixMax12d EdgeEdgeCollision::mollifier_hessian(
    const VectorMax12d& positions, const double _eps_x) const
{
    assert(positions.size() == 12);
    return edge_edge_mollifier_hessian(
        positions.segment<3>(0), positions.segment<3>(3),
        positions.segment<3>(6), positions.segment<3>(9), _eps_x);
}

Vector12d EdgeEdgeCollision::mollifier_gradient_wrt_x(
    const VectorMax12d& rest_positions, const VectorMax12d& positions) const
{
    assert(rest_positions.size() == 12);
    assert(positions.size() == 12);
    return edge_edge_mollifier_gradient_wrt_x(
        rest_positions.segment<3>(0), rest_positions.segment<3>(3),
        rest_positions.segment<3>(6), rest_positions.segment<3>(9),
        positions.segment<3>(0), positions.segment<3>(3),
        positions.segment<3>(6), positions.segment<3>(9));
}

Matrix12d EdgeEdgeCollision::mollifier_gradient_jacobian_wrt_x(
    const VectorMax12d& rest_positions, const VectorMax12d& positions) const
{
    assert(rest_positions.size() == 12);
    assert(positions.size() == 12);
    return edge_edge_mollifier_gradient_jacobian_wrt_x(
        rest_positions.segment<3>(0), rest_positions.segment<3>(3),
        rest_positions.segment<3>(6), rest_positions.segment<3>(9),
        positions.segment<3>(0), positions.segment<3>(3),
        positions.segment<3>(6), positions.segment<3>(9));
}

bool EdgeEdgeCollision::operator==(const EdgeEdgeCollision& other) const
{
    return EdgeEdgeCandidate::operator==(other) && dtype == other.dtype;
}

bool EdgeEdgeCollision::operator!=(const EdgeEdgeCollision& other) const
{
    return !(*this == other);
}

bool EdgeEdgeCollision::operator<(const EdgeEdgeCollision& other) const
{
    if (EdgeEdgeCandidate::operator==(other)) {
        return dtype < other.dtype;
    }
    return EdgeEdgeCandidate::operator<(other);
}

} // namespace ipc
