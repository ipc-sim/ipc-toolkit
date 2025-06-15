#include "edge2.hpp"

namespace ipc {
Edge2::Edge2(
    const long& id,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const VectorMax3d& d,
    const ParameterType& param)
    : Primitive(id, param)
{
    _vert_ids = { { mesh.edges()(id, 0), mesh.edges()(id, 1) } };

    is_active_ = (mesh.is_orient_vertex(_vert_ids[0]) && mesh.is_orient_vertex(_vert_ids[1]))
        || Math<double>::cross2(
               d, vertices.row(_vert_ids[1]) - vertices.row(_vert_ids[0]))
            > 0;
}

int Edge2::n_vertices() const { return n_edge_neighbors_2d; }

double Edge2::potential(const Vector2d& d, const Vector4d& x) const
{
    return (x.tail<2>() - x.head<2>()).norm();
}

Vector6d Edge2::grad(const Vector2d& d, const Vector4d& x) const
{
    const Vector2d t = x.tail<2>() - x.head<2>();
    const double len = t.norm();
    Vector6d g;
    g.setZero();
    g.segment<2>(2) = -t / len;
    g.segment<2>(4) = t / len;
    return g;
}

Matrix6d Edge2::hessian(const Vector2d& d, const Vector4d& x) const
{
    Matrix6d h;
    h.setZero();
#ifdef DERIVATIVES_WITH_AUTODIFF
    DiffScalarBase::setVariableCount(4);
    using T = ADHessian<4>;
    auto xAD = slice_positions<T, 2, 2>(x);
    h.block<4, 4>(2, 2) = (xAD.row(0) - xAD.row(1)).norm().getHessian();
#else
    const Vector2d t = x.tail<2>() - x.head<2>();
    const double norm = t.norm();
    h.block<2, 2>(2, 2) =
        (Matrix2d::Identity() - t * (1. / norm / norm) * t.transpose()) / norm;
    h.block<2, 2>(4, 4) = h.block<2, 2>(2, 2);
    h.block<2, 2>(2, 4) = -h.block<2, 2>(2, 2);
    h.block<2, 2>(4, 2) = -h.block<2, 2>(2, 2);
#endif
    return h;
}
} // namespace ipc
