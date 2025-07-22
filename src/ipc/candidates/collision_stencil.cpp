#include "collision_stencil.hpp"

namespace ipc {

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

    VectorMax3d n = compute_unnormalized_normal(positions);

    MatrixMax<double, 3, 12> dn =
        compute_unnormalized_normal_jacobian(positions);

#if true
    // Derivative of normalization (n̂ = n / ‖n‖)
    const double n_norm = n.norm();
    n /= n_norm; // n̂

    MatrixMax3d A =
        (MatrixMax3d::Identity(dim, dim) - n * n.transpose()) / n_norm;

    dn = A * dn;
#endif

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
