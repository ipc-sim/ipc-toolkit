#include "edge.hpp"

#include <ipc/config.hpp>

namespace ipc {

template <int DIM>
Edge<DIM>::Edge(
    const index_t id,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const VectorMax3d& d,
    const SmoothContactParameters& params)
    : Primitive(id, params)
{
    m_vertex_ids = { { mesh.edges()(id, 0), mesh.edges()(id, 1) } };

    m_is_active = (mesh.is_orient_vertex(m_vertex_ids[0])
                   && mesh.is_orient_vertex(m_vertex_ids[1]))
        || Math<double>::cross2(
               d, vertices.row(m_vertex_ids[1]) - vertices.row(m_vertex_ids[0]))
            > 0;
}

template <> int Edge<2>::n_vertices() const { return N_EDGE_NEIGHBORS_2D; }
template <> int Edge<3>::n_vertices() const { return N_EDGE_NEIGHBORS_3D; }

template <>
double
Edge<2>::potential(Eigen::ConstRef<DVector> d, Eigen::ConstRef<XVector> x) const
{
    return (x.tail<2>() - x.head<2>()).norm();
}

template <>
Vector6d
Edge<2>::grad(Eigen::ConstRef<DVector> d, Eigen::ConstRef<XVector> x) const
{
    const Eigen::Vector2d t = x.tail<2>() - x.head<2>();
    const double len = t.norm();
    Vector6d g;
    g.setZero();
    g.segment<2>(2) = -t / len;
    g.segment<2>(4) = t / len;
    return g;
}

template <>
Matrix6d
Edge<2>::hessian(Eigen::ConstRef<DVector> d, Eigen::ConstRef<XVector> x) const
{
    Matrix6d h;
    h.setZero();
#ifdef IPC_TOOLKIT_DEBUG_AUTODIFF
    ScalarBase::setVariableCount(4);
    using T = ADHessian<4>;
    auto xAD = slice_positions<T, 2, 2>(x);
    h.block<4, 4>(2, 2) = (xAD.row(0) - xAD.row(1)).norm().Hess;
#else
    const Eigen::Vector2d t = x.tail<2>() - x.head<2>();
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
