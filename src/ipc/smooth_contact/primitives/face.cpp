#include "face.hpp"

#include <ipc/utils/autodiff_types.hpp>

namespace ipc {
// namespace {
//     template <typename T>
//     std::array<Vector3<T>, 3> double_to_autodiff(
//         Eigen::ConstRef<Eigen::Vector3d> v0,
//         Eigen::ConstRef<Eigen::Vector3d> v1,
//         Eigen::ConstRef<Eigen::Vector3d> v2)
//     {
//         Vector3<T> v0_, v1_, v2_;
//         for (int d = 0; d < 3; d++)
//         {
//             v0_(d) = T(  d, v0(d));
//             v1_(d) = T(3+d, v1(d));
//             v2_(d) = T(6+d, v2(d));
//         }
//         return {{v0_, v1_, v2_}};
//     }
// }
Face::Face(
    const index_t id,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const VectorMax3d& d,
    const SmoothContactParameters& params)
    : Primitive(id, params)
{
    m_vertex_ids = { { mesh.faces()(id, 0), mesh.faces()(id, 1),
                       mesh.faces()(id, 2) } };
    Vector3d a = vertices.row(m_vertex_ids[1]) - vertices.row(m_vertex_ids[0]);
    Vector3d b = vertices.row(m_vertex_ids[2]) - vertices.row(m_vertex_ids[0]);

    bool orientable = mesh.is_orient_vertex(m_vertex_ids[0])
        && mesh.is_orient_vertex(m_vertex_ids[1])
        && mesh.is_orient_vertex(m_vertex_ids[2]);
    m_is_active = !orientable || a.cross(b).dot(d) > 0;
}
int Face::n_vertices() const { return N_FACE_NEIGHBORS_3D; }
double Face::potential(const Vector3d& d, const Vector9d& x) const
{
    return smooth_face_term<double>(x.head<3>(), x.segment<3>(3), x.tail<3>());
}
Vector12d Face::grad(const Vector3d& d, const Vector9d& x) const
{
    Vector12d g;
    g.setZero();
    auto X = slice_positions<TinyADGrad<9>, 3, 3>(x);
    g.tail<9>() =
        smooth_face_term<TinyADGrad<9>>(X.row(0), X.row(1), X.row(2)).grad;
    return g;
}
Matrix12d Face::hessian(const Vector3d& d, const Vector9d& x) const
{
    Matrix12d h;
    h.setZero();
    auto X = slice_positions<TinyADHessian<9>, 3, 3>(x);
    h.bottomRightCorner<9, 9>() =
        smooth_face_term<TinyADHessian<9>>(X.row(0), X.row(1), X.row(2))
            .Hess;
    return h;
}
} // namespace ipc
