#include "collision_stencil.hpp"

#include <ipc/geometry/normal.hpp>

namespace ipc {

VectorMax3d CollisionStencil::compute_distance_vector(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    VectorMax4d _; // Unused output
    return compute_distance_vector(positions, _);
}

VectorMax3d CollisionStencil::compute_distance_vector(
    Eigen::ConstRef<VectorMax12d> positions, VectorMax4d& coeffs) const
{
    const int n = num_vertices();
    const int d = dim(positions.size());
    coeffs = compute_coefficients(positions);
    VectorMax3d dv = VectorMax3d::Zero(d);
    for (int i = 0; i < n; i++) {
        dv += coeffs[i] * positions.segment(d * i, d);
    }
    return dv;
}

MatrixMax<double, 12, 3> CollisionStencil::compute_distance_vector_jacobian(
    Eigen::ConstRef<VectorMax12d> positions) const
{
    const int n = num_vertices();
    const int d = dim(positions.size());
    const int ndof = n * d;
    const VectorMax4d c = compute_coefficients(positions);
    // ∂t/∂x is ndof × dim
    // Each block row i is cᵢ * I_{d×d}
    MatrixMax<double, 12, 3> J;
    J.setZero(ndof, d);
    for (int i = 0; i < n; i++) {
        J.block(i * d, 0, d, d).diagonal().array() = c[i];
    }
    return J;
}

VectorMax12d CollisionStencil::diag_distance_vector_outer(
    Eigen::ConstRef<VectorMax4d> coeffs, const int d)
{
    const int n = coeffs.size();
    VectorMax12d diag(n * d);
    for (int i = 0; i < n; i++) {
        const double ci2 = coeffs[i] * coeffs[i];
        diag.segment(i * d, d).setConstant(ci2);
    }
    return diag;
}

VectorMax12d CollisionStencil::diag_distance_vector_t_outer(
    Eigen::ConstRef<VectorMax4d> coeffs,
    Eigen::ConstRef<VectorMax3d> distance_vector)
{
    const int n = coeffs.size();
    const int d = distance_vector.size();
    const VectorMax3d t2 = distance_vector.array().square();
    VectorMax12d diag(n * d);
    for (int i = 0; i < n; i++) {
        diag.segment(i * d, d).array() = (coeffs[i] * coeffs[i]) * t2;
    }
    return diag;
}

VectorMax3d CollisionStencil::contract_distance_vector_jacobian(
    Eigen::ConstRef<VectorMax4d> coeffs,
    Eigen::ConstRef<VectorMax12d> p,
    const int d)
{
    const int n = coeffs.size();
    VectorMax3d result = VectorMax3d::Zero(d);
    for (int i = 0; i < n; i++) {
        result += coeffs[i] * p.segment(d * i, d);
    }
    return result;
}

VectorMax3d CollisionStencil::compute_normal(
    Eigen::ConstRef<VectorMax12d> positions,
    bool flip_if_negative,
    double* sign) const
{
    const int dim = this->dim(positions.size());

    VectorMax3d n = compute_unnormalized_normal(positions).normalized();

    if (sign != nullptr) {
        *sign = (positions.head(dim) - positions.tail(dim)).dot(n) < 0 ? -1 : 1;
    }

    // Flip the normal if the point is on the negative side.
    // Any point on the second object will do, so we use the last point.
    if (flip_if_negative
        && (positions.head(dim) - positions.tail(dim)).dot(n) < 0) {
        n *= -1;
    }

    return n;
}

MatrixMax<double, 3, 12> CollisionStencil::compute_normal_jacobian(
    Eigen::ConstRef<VectorMax12d> positions, bool flip_if_negative) const
{
    const int dim = this->dim(positions.size());

    const VectorMax3d n = compute_unnormalized_normal(positions);

    MatrixMax<double, 3, 12> dn = normalization_jacobian(n)
        * compute_unnormalized_normal_jacobian(positions);

    if (flip_if_negative
        && (positions.head(dim) - positions.tail(dim)).dot(n) < 0) {
        dn *= -1;
    }

    return dn;
}

std::ostream& CollisionStencil::write_ccd_query(
    std::ostream& out,
    Eigen::ConstRef<VectorMax12d> vertices_t0,
    Eigen::ConstRef<VectorMax12d> vertices_t1) const
{
    assert(vertices_t0.size() == vertices_t1.size());

    const int dim = vertices_t0.size() / num_vertices();
    assert(vertices_t0.size() % num_vertices() == 0);

    for (int i = 0; i < num_vertices(); i++) {
        out << vertices_t0.segment(dim * i, dim)
                   .transpose()
                   .format(OBJ_VERTEX_FORMAT);
    }

    for (int i = 0; i < num_vertices(); i++) {
        out << vertices_t1.segment(dim * i, dim)
                   .transpose()
                   .format(OBJ_VERTEX_FORMAT);
    }

    return out;
}

} // namespace ipc
