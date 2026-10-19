#include "esp_collision_template.hpp"

#include <ipc/distance/distance_type_exact.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <algorithm>

namespace ipc {

// ---- type ----

template <>
ESPCollisionType ESPCollisionTemplate<Vertex3, Vertex3>::type() const
{
    return ESPCollisionType::VERTEX_VERTEX;
}
template <>
ESPCollisionType ESPCollisionTemplate<Edge3P1, Vertex3>::type() const
{
    return ESPCollisionType::EDGE_VERTEX;
}
template <>
ESPCollisionType ESPCollisionTemplate<Face3P1, Vertex3>::type() const
{
    return ESPCollisionType::FACE_VERTEX;
}
template <>
ESPCollisionType ESPCollisionTemplate<Vertex2, Vertex2>::type() const
{
    return ESPCollisionType::VERTEX_VERTEX;
}
template <>
ESPCollisionType ESPCollisionTemplate<Vertex2, Edge2P1>::type() const
{
    return ESPCollisionType::EDGE_VERTEX;
}

// ---- name ----

template <> std::string ESPCollisionTemplate<Vertex3, Vertex3>::name() const
{
    return "vv_3d";
}
template <> std::string ESPCollisionTemplate<Edge3P1, Vertex3>::name() const
{
    return "ev_3d";
}
template <> std::string ESPCollisionTemplate<Face3P1, Vertex3>::name() const
{
    return "fv_3d";
}
template <> std::string ESPCollisionTemplate<Vertex2, Vertex2>::name() const
{
    return "vv_2d_pt";
}
template <> std::string ESPCollisionTemplate<Vertex2, Edge2P1>::name() const
{
    return "ev_2d_pt";
}

// ---- constructors ----

template <typename PrimitiveA, typename PrimitiveB>
ESPCollisionTemplate<PrimitiveA, PrimitiveB>::ESPCollisionTemplate(
    index_t _primitive0, index_t _primitive1, const CollisionMesh& mesh)
    : m_primitive_a(_primitive0, mesh)
    , m_primitive_b(_primitive1, mesh)
{
}

template <>
ESPCollisionTemplate<Vertex3, Vertex3>::ESPCollisionTemplate(
    index_t _primitive0, index_t _primitive1, const CollisionMesh& mesh)
    : m_primitive_a(std::min(_primitive0, _primitive1), mesh)
    , m_primitive_b(std::max(_primitive0, _primitive1), mesh)
{
}

// ---- vertex_id ----

template <typename PrimitiveA, typename PrimitiveB>
index_t ESPCollisionTemplate<PrimitiveA, PrimitiveB>::vertex_id(index_t i) const
{
    if (i < static_cast<index_t>(m_primitive_a.n_vertices())) {
        return m_primitive_a.vertex_id(i);
    }
    i -= m_primitive_a.n_vertices();
    assert(static_cast<index_t>(m_primitive_b.n_vertices()) > i);
    return m_primitive_b.vertex_id(i);
}

// ---- generic stubs ----

template <typename PrimitiveA, typename PrimitiveB>
double ESPCollisionTemplate<PrimitiveA, PrimitiveB>::operator()(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> /*positions*/,
    const ESPParameters& /*params*/) const
{
    return 0;
}

template <typename PrimitiveA, typename PrimitiveB>
auto ESPCollisionTemplate<PrimitiveA, PrimitiveB>::gradient(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> /*positions*/,
    const ESPParameters& /*params*/) const -> VectorMax<double, ELEMENT_SIZE>
{
    return VectorMax<double, ELEMENT_SIZE>::Zero(n_dofs());
}

template <typename PrimitiveA, typename PrimitiveB>
auto ESPCollisionTemplate<PrimitiveA, PrimitiveB>::hessian(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> /*positions*/,
    const ESPParameters& /*params*/) const
    -> MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>
{
    return MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>::Zero(
        n_dofs(), n_dofs());
}

template <typename PrimitiveA, typename PrimitiveB>
double ESPCollisionTemplate<PrimitiveA, PrimitiveB>::compute_distance(
    Eigen::ConstRef<Eigen::MatrixXd> /*vertices*/) const
{
    log_and_throw_error("Not implemented");
    return 0;
}

template <typename PrimitiveA, typename PrimitiveB>
std::pair<double, double>
ESPCollisionTemplate<PrimitiveA, PrimitiveB>::operator_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> /*positions*/,
    const ESPParameters& /*params*/,
    const NearFarBarrier* /*nf_barrier*/) const
{
    return { 0.0, 0.0 };
}

template <typename PrimitiveA, typename PrimitiveB>
auto ESPCollisionTemplate<PrimitiveA, PrimitiveB>::gradient_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& /*params*/,
    const NearFarBarrier* /*nf_barrier*/) const -> std::
    pair<VectorMax<double, ELEMENT_SIZE>, VectorMax<double, ELEMENT_SIZE>>
{
    const auto n = positions.size();
    return { VectorMax<double, ELEMENT_SIZE>::Zero(n),
             VectorMax<double, ELEMENT_SIZE>::Zero(n) };
}

template <typename PrimitiveA, typename PrimitiveB>
auto ESPCollisionTemplate<PrimitiveA, PrimitiveB>::hessian_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& /*params*/,
    const NearFarBarrier* /*nf_barrier*/) const
    -> std::pair<
        MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>,
        MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>>
{
    const auto n = positions.size();
    return { MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>::Zero(n, n),
             MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>::Zero(n, n) };
}

// ---- 3D specializations ----

template <>
double ESPCollisionTemplate<Vertex3, Vertex3>::compute_distance(
    Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    const int n_verts = vertices.rows();
    if (n_verts > vertex_id(0) && n_verts > vertex_id(1)) {
        return point_point_distance(
            vertices.row(vertex_id(0)),
            vertices.row(vertex_id(n_vertices_a())));
    }
    return std::numeric_limits<double>::max();
}

template <>
double ESPCollisionTemplate<Edge3P1, Vertex3>::compute_distance(
    Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    const int n_verts = vertices.rows();
    if (n_verts > vertex_id(0) && n_verts > vertex_id(1)
        && n_verts > vertex_id(2)) {
        return point_edge_distance(
            vertices.row(vertex_id(n_vertices_a())), vertices.row(vertex_id(0)),
            vertices.row(vertex_id(1)));
    }
    return std::numeric_limits<double>::max();
}

template <>
double ESPCollisionTemplate<Face3P1, Vertex3>::compute_distance(
    Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    const int n_verts = vertices.rows();
    if (n_verts > vertex_id(0) && n_verts > vertex_id(1)
        && n_verts > vertex_id(2) && n_verts > vertex_id(3)) {
        return point_triangle_distance(
            vertices.row(vertex_id(3)), vertices.row(vertex_id(0)),
            vertices.row(vertex_id(1)), vertices.row(vertex_id(2)));
    }
    return std::numeric_limits<double>::max();
}

template <>
double ESPCollisionTemplate<Vertex3, Vertex3>::operator()(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params) const
{
    const double dist =
        (positions.template head<3>() - positions.template segment<3>(3))
            .norm();
    const double eps = params.dhat;
    return (*params.barrier)(dist, eps);
}

template <>
double ESPCollisionTemplate<Edge3P1, Vertex3>::operator()(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params) const
{
    assert(
        point_edge_distance_type_exact(
            positions.template segment<3>(6), positions.template head<3>(),
            positions.template segment<3>(3))
        == PointEdgeDistanceType::P_E);
    const double eps = params.dhat;
    // Edge3P1-Vertex3 is constructed only at interior P_E (see
    // ESPCollisionsBuilder<3>::reduce_point_edge_collision).
    const double dist = sqrt(point_edge_distance(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3)));
    return (*params.barrier)(dist, eps);
}

template <>
double ESPCollisionTemplate<Face3P1, Vertex3>::operator()(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params) const
{
    assert(
        point_triangle_distance_type_exact(
            positions.template segment<3>(9), positions.template head<3>(),
            positions.template segment<3>(3), positions.template segment<3>(6))
        == PointTriangleDistanceType::P_T);
    const double eps = params.dhat;
    // Face3P1-Vertex3 is constructed only at interior P_T (see
    // ESPCollisionsBuilder<3>::reduce_point_triangle_collision).
    const double dist = sqrt(point_triangle_distance(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6)));
    return (*params.barrier)(dist, eps);
}

template <>
auto ESPCollisionTemplate<Vertex3, Vertex3>::gradient(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params) const -> VectorMax<double, ELEMENT_SIZE>
{
    assert(positions.size() == 6);
    const double dist =
        (positions.template head<3>() - positions.template tail<3>()).norm();
    const double eps = params.dhat;
    const double deriv =
        params.barrier->first_derivative(dist, eps) / (dist * 2.);
    Vector6d grad =
        deriv
        * point_point_distance_gradient(
            positions.template head<3>(), positions.template tail<3>());
    return grad;
}

template <>
auto ESPCollisionTemplate<Edge3P1, Vertex3>::gradient(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params) const -> VectorMax<double, ELEMENT_SIZE>
{
    assert(positions.size() == 9);
    assert(
        point_edge_distance_type_exact(
            positions.template segment<3>(6), positions.template head<3>(),
            positions.template segment<3>(3))
        == PointEdgeDistanceType::P_E);
    // Edge3P1-Vertex3 is constructed only at interior P_E.
    constexpr auto dtype = PointEdgeDistanceType::P_E;
    const double dist = sqrt(point_edge_distance(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3), dtype));
    const double eps = params.dhat;
    const double deriv =
        params.barrier->first_derivative(dist, eps) / (dist * 2.);
    Vector9d grad = point_edge_distance_gradient(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3), dtype);
    grad *= deriv;
    grad = grad({ 3, 4, 5, 6, 7, 8, 0, 1, 2 }).eval();
    return grad;
}

template <>
auto ESPCollisionTemplate<Face3P1, Vertex3>::gradient(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params) const -> VectorMax<double, ELEMENT_SIZE>
{
    assert(positions.size() == 12);
    assert(
        point_triangle_distance_type_exact(
            positions.template segment<3>(9), positions.template head<3>(),
            positions.template segment<3>(3), positions.template segment<3>(6))
        == PointTriangleDistanceType::P_T);
    // Face3P1-Vertex3 is constructed only at interior P_T.
    constexpr auto dtype = PointTriangleDistanceType::P_T;
    const double dist = sqrt(point_triangle_distance(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6),
        dtype));
    const double eps = params.dhat;
    const double deriv =
        params.barrier->first_derivative(dist, eps) / (dist * 2.);
    Vector12d grad = point_triangle_distance_gradient(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6),
        dtype);
    grad *= deriv;
    grad = grad({ 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2 }).eval();
    return grad;
}

template <>
auto ESPCollisionTemplate<Vertex3, Vertex3>::hessian(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params) const
    -> MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>
{
    assert(positions.size() == 6);
    const double dist =
        (positions.template head<3>() - positions.template tail<3>()).norm();
    const double eps = params.dhat;
    double deriv1 = params.barrier->first_derivative(dist, eps);
    double deriv2 = params.barrier->second_derivative(dist, eps);
    deriv2 = deriv2 / (4 * dist * dist) - deriv1 / (4 * dist * dist * dist);
    deriv1 /= (2 * dist);
    const Vector6d g = point_point_distance_gradient(
        positions.template head<3>(), positions.template tail<3>());
    const Matrix6d h = point_point_distance_hessian(
        positions.template head<3>(), positions.template tail<3>());
    return g * deriv2 * g.transpose() + h * deriv1;
}

template <>
auto ESPCollisionTemplate<Edge3P1, Vertex3>::hessian(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params) const
    -> MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>
{
    assert(positions.size() == 9);
    assert(
        point_edge_distance_type_exact(
            positions.template segment<3>(6), positions.template head<3>(),
            positions.template segment<3>(3))
        == PointEdgeDistanceType::P_E);
    // Edge3P1-Vertex3 is constructed only at interior P_E.
    constexpr auto dtype = PointEdgeDistanceType::P_E;
    const double dist = sqrt(point_edge_distance(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3), dtype));
    const double eps = params.dhat;
    double deriv1 = params.barrier->first_derivative(dist, eps);
    double deriv2 = params.barrier->second_derivative(dist, eps);
    deriv2 = deriv2 / (4 * dist * dist) - deriv1 / (4 * dist * dist * dist);
    deriv1 /= (2 * dist);
    const Vector9d g = point_edge_distance_gradient(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3), dtype);
    const Matrix9d h = point_edge_distance_hessian(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3), dtype);
    Matrix9d hess = g * deriv2 * g.transpose() + h * deriv1;
    std::vector<int> reorder { 3, 4, 5, 6, 7, 8, 0, 1, 2 };
    return hess(reorder, reorder);
}

template <>
auto ESPCollisionTemplate<Face3P1, Vertex3>::hessian(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params) const
    -> MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>
{
    assert(positions.size() == 12);
    assert(
        point_triangle_distance_type_exact(
            positions.template segment<3>(9), positions.template head<3>(),
            positions.template segment<3>(3), positions.template segment<3>(6))
        == PointTriangleDistanceType::P_T);
    // Face3P1-Vertex3 is constructed only at interior P_T.
    constexpr auto dtype = PointTriangleDistanceType::P_T;
    const double dist = sqrt(point_triangle_distance(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6),
        dtype));
    const double eps = params.dhat;
    double deriv1 = params.barrier->first_derivative(dist, eps);
    double deriv2 = params.barrier->second_derivative(dist, eps);
    deriv2 = deriv2 / (4 * dist * dist) - deriv1 / (4 * dist * dist * dist);
    deriv1 /= (2 * dist);
    const Vector12d g = point_triangle_distance_gradient(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6),
        dtype);
    const Matrix12d h = point_triangle_distance_hessian(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6),
        dtype);
    Matrix12d hess = g * deriv2 * g.transpose() + h * deriv1;
    std::vector<int> reorder { 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2 };
    return hess(reorder, reorder);
}

// ---- NearFarBarrier specializations (3D only) ----

template <>
std::pair<double, double>
ESPCollisionTemplate<Vertex3, Vertex3>::operator_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const NearFarBarrier* nf_barrier) const
{
    const double dist =
        (positions.template head<3>() - positions.template segment<3>(3))
            .norm();
    const double eps = params.dhat;
    return { nf_barrier->near_value(dist, eps),
             nf_barrier->far_value(dist, eps) };
}

template <>
std::pair<double, double>
ESPCollisionTemplate<Edge3P1, Vertex3>::operator_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const NearFarBarrier* nf_barrier) const
{
    const double dist = sqrt(point_edge_distance(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3)));
    const double eps = params.dhat;
    return { nf_barrier->near_value(dist, eps),
             nf_barrier->far_value(dist, eps) };
}

template <>
std::pair<double, double>
ESPCollisionTemplate<Face3P1, Vertex3>::operator_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const NearFarBarrier* nf_barrier) const
{
    const double dist = sqrt(point_triangle_distance(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6)));
    const double eps = params.dhat;
    return { nf_barrier->near_value(dist, eps),
             nf_barrier->far_value(dist, eps) };
}

template <>
std::pair<
    VectorMax<double, ESPCollision::ELEMENT_SIZE>,
    VectorMax<double, ESPCollision::ELEMENT_SIZE>>
ESPCollisionTemplate<Vertex3, Vertex3>::gradient_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const NearFarBarrier* nf_barrier) const
{
    assert(positions.size() == 6);
    const double dist =
        (positions.template head<3>() - positions.template tail<3>()).norm();
    const double eps = params.dhat;
    const double deriv_near =
        nf_barrier->first_derivative_near(dist, eps) / (dist * 2.);
    const double deriv_far =
        nf_barrier->first_derivative_far(dist, eps) / (dist * 2.);
    Vector6d g = point_point_distance_gradient(
        positions.template head<3>(), positions.template tail<3>());
    return { deriv_near * g, deriv_far * g };
}

template <>
std::pair<
    VectorMax<double, ESPCollision::ELEMENT_SIZE>,
    VectorMax<double, ESPCollision::ELEMENT_SIZE>>
ESPCollisionTemplate<Edge3P1, Vertex3>::gradient_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const NearFarBarrier* nf_barrier) const
{
    assert(positions.size() == 9);
    auto dtype = point_edge_distance_type_exact(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3));
    const double dist = sqrt(point_edge_distance(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3), dtype));
    const double eps = params.dhat;
    const double deriv_near =
        nf_barrier->first_derivative_near(dist, eps) / (dist * 2.);
    const double deriv_far =
        nf_barrier->first_derivative_far(dist, eps) / (dist * 2.);
    Vector9d g = point_edge_distance_gradient(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3), dtype);
    g = g({ 3, 4, 5, 6, 7, 8, 0, 1, 2 }).eval();
    return { deriv_near * g, deriv_far * g };
}

template <>
std::pair<
    VectorMax<double, ESPCollision::ELEMENT_SIZE>,
    VectorMax<double, ESPCollision::ELEMENT_SIZE>>
ESPCollisionTemplate<Face3P1, Vertex3>::gradient_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const NearFarBarrier* nf_barrier) const
{
    assert(positions.size() == 12);
    auto dtype = point_triangle_distance_type_exact(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6));
    const double dist = sqrt(point_triangle_distance(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6),
        dtype));
    const double eps = params.dhat;
    const double deriv_near =
        nf_barrier->first_derivative_near(dist, eps) / (dist * 2.);
    const double deriv_far =
        nf_barrier->first_derivative_far(dist, eps) / (dist * 2.);
    Vector12d g = point_triangle_distance_gradient(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6),
        dtype);
    g = g({ 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2 }).eval();
    return { deriv_near * g, deriv_far * g };
}

template <>
std::pair<
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>,
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>>
ESPCollisionTemplate<Vertex3, Vertex3>::hessian_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const NearFarBarrier* nf_barrier) const
{
    assert(positions.size() == 6);
    const double dist =
        (positions.template head<3>() - positions.template tail<3>()).norm();
    const double eps = params.dhat;
    double deriv1_near = nf_barrier->first_derivative_near(dist, eps);
    double deriv2_near = nf_barrier->second_derivative_near(dist, eps);
    double deriv1_far = nf_barrier->first_derivative_far(dist, eps);
    double deriv2_far = nf_barrier->second_derivative_far(dist, eps);
    deriv2_near = deriv2_near / (4 * dist * dist)
        - deriv1_near / (4 * dist * dist * dist);
    deriv2_far =
        deriv2_far / (4 * dist * dist) - deriv1_far / (4 * dist * dist * dist);
    deriv1_near /= (2 * dist);
    deriv1_far /= (2 * dist);
    const Vector6d g = point_point_distance_gradient(
        positions.template head<3>(), positions.template tail<3>());
    const Matrix6d h = point_point_distance_hessian(
        positions.template head<3>(), positions.template tail<3>());
    return { g * deriv2_near * g.transpose() + h * deriv1_near,
             g * deriv2_far * g.transpose() + h * deriv1_far };
}

template <>
std::pair<
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>,
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>>
ESPCollisionTemplate<Edge3P1, Vertex3>::hessian_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const NearFarBarrier* nf_barrier) const
{
    assert(positions.size() == 9);
    auto dtype = point_edge_distance_type_exact(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3));
    const double dist = sqrt(point_edge_distance(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3), dtype));
    const double eps = params.dhat;
    double deriv1_near = nf_barrier->first_derivative_near(dist, eps);
    double deriv2_near = nf_barrier->second_derivative_near(dist, eps);
    double deriv1_far = nf_barrier->first_derivative_far(dist, eps);
    double deriv2_far = nf_barrier->second_derivative_far(dist, eps);
    deriv2_near = deriv2_near / (4 * dist * dist)
        - deriv1_near / (4 * dist * dist * dist);
    deriv2_far =
        deriv2_far / (4 * dist * dist) - deriv1_far / (4 * dist * dist * dist);
    deriv1_near /= (2 * dist);
    deriv1_far /= (2 * dist);
    const Vector9d g = point_edge_distance_gradient(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3), dtype);
    const Matrix9d h = point_edge_distance_hessian(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3), dtype);
    const Matrix9d hess_near =
        g * deriv2_near * g.transpose() + h * deriv1_near;
    const Matrix9d hess_far = g * deriv2_far * g.transpose() + h * deriv1_far;
    const std::array<int, 9> reorder { { 3, 4, 5, 6, 7, 8, 0, 1, 2 } };
    return { hess_near(reorder, reorder), hess_far(reorder, reorder) };
}

template <>
std::pair<
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>,
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>>
ESPCollisionTemplate<Face3P1, Vertex3>::hessian_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const NearFarBarrier* nf_barrier) const
{
    assert(positions.size() == 12);
    auto dtype = point_triangle_distance_type_exact(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6));
    const double dist = sqrt(point_triangle_distance(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6),
        dtype));
    const double eps = params.dhat;
    double deriv1_near = nf_barrier->first_derivative_near(dist, eps);
    double deriv2_near = nf_barrier->second_derivative_near(dist, eps);
    double deriv1_far = nf_barrier->first_derivative_far(dist, eps);
    double deriv2_far = nf_barrier->second_derivative_far(dist, eps);
    deriv2_near = deriv2_near / (4 * dist * dist)
        - deriv1_near / (4 * dist * dist * dist);
    deriv2_far =
        deriv2_far / (4 * dist * dist) - deriv1_far / (4 * dist * dist * dist);
    deriv1_near /= (2 * dist);
    deriv1_far /= (2 * dist);
    const Vector12d g = point_triangle_distance_gradient(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6),
        dtype);
    const Matrix12d h = point_triangle_distance_hessian(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6),
        dtype);
    const Matrix12d hess_near =
        g * deriv2_near * g.transpose() + h * deriv1_near;
    const Matrix12d hess_far = g * deriv2_far * g.transpose() + h * deriv1_far;
    const std::array<int, 12> reorder { { 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1,
                                          2 } };
    return { hess_near(reorder, reorder), hess_far(reorder, reorder) };
}

// ---- 2D specializations ----
// positions layout VV: [q_x, q_y, v_x, v_y]
// positions layout VE: [q_x, q_y, e0_x, e0_y, e1_x, e1_y]

template <>
double ESPCollisionTemplate<Vertex2, Vertex2>::compute_distance(
    Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    const int n = vertices.rows();
    if (vertex_id(0) >= n || vertex_id(1) >= n) {
        return std::numeric_limits<double>::max();
    }
    return point_point_distance(
        vertices.row(vertex_id(0)), vertices.row(vertex_id(1)));
}

template <>
double ESPCollisionTemplate<Vertex2, Edge2P1>::compute_distance(
    Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    const int n = vertices.rows();
    if (vertex_id(0) >= n || vertex_id(1) >= n || vertex_id(2) >= n) {
        return std::numeric_limits<double>::max();
    }
    return point_edge_distance(
        vertices.row(vertex_id(0)), vertices.row(vertex_id(1)),
        vertices.row(vertex_id(2)));
}

template <>
double ESPCollisionTemplate<Vertex2, Vertex2>::operator()(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params) const
{
    const double dist =
        (positions.template head<2>() - positions.template tail<2>()).norm();
    const double eps = params.dhat;
    return (*params.barrier)(dist, eps);
}

template <>
double ESPCollisionTemplate<Vertex2, Edge2P1>::operator()(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params) const
{
    assert(
        point_edge_distance_type_exact(
            positions.template head<2>(), positions.template segment<2>(2),
            positions.template segment<2>(4))
        == PointEdgeDistanceType::P_E);
    const double eps = params.dhat;
    // Vertex2-Edge2P1 is constructed only at interior P_E (the 2D edge-QP
    // builder routes endpoint cases to Vertex2-Vertex2).
    const double dist = std::sqrt(point_edge_distance(
        positions.template head<2>(), positions.template segment<2>(2),
        positions.template segment<2>(4)));
    return (*params.barrier)(dist, eps);
}

template <>
auto ESPCollisionTemplate<Vertex2, Vertex2>::gradient(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params) const -> VectorMax<double, ELEMENT_SIZE>
{
    const double dist =
        (positions.template head<2>() - positions.template tail<2>()).norm();
    const double eps = params.dhat;
    const double deriv =
        params.barrier->first_derivative(dist, eps) / (dist * 2.0);
    const VectorMax6d g = point_point_distance_gradient(
        positions.template head<2>(), positions.template tail<2>());
    return deriv * g;
}

template <>
auto ESPCollisionTemplate<Vertex2, Edge2P1>::gradient(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params) const -> VectorMax<double, ELEMENT_SIZE>
{
    assert(
        point_edge_distance_type_exact(
            positions.template head<2>(), positions.template segment<2>(2),
            positions.template segment<2>(4))
        == PointEdgeDistanceType::P_E);
    // Vertex2-Edge2P1 is constructed only at interior P_E.
    constexpr auto dtype = PointEdgeDistanceType::P_E;
    const double dist = std::sqrt(point_edge_distance(
        positions.template head<2>(), positions.template segment<2>(2),
        positions.template segment<2>(4), dtype));
    const double eps = params.dhat;
    const double deriv =
        params.barrier->first_derivative(dist, eps) / (dist * 2.0);
    const VectorMax9d g = point_edge_distance_gradient(
        positions.template head<2>(), positions.template segment<2>(2),
        positions.template segment<2>(4), dtype);
    return deriv * g;
}

template <>
auto ESPCollisionTemplate<Vertex2, Vertex2>::hessian(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params) const
    -> MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>
{
    const double dist =
        (positions.template head<2>() - positions.template tail<2>()).norm();
    const double eps = params.dhat;
    double deriv1 = params.barrier->first_derivative(dist, eps);
    double deriv2 = params.barrier->second_derivative(dist, eps);
    deriv2 = deriv2 / (4 * dist * dist) - deriv1 / (4 * dist * dist * dist);
    deriv1 /= (2 * dist);
    const VectorMax6d g = point_point_distance_gradient(
        positions.template head<2>(), positions.template tail<2>());
    const MatrixMax6d H = point_point_distance_hessian(
        positions.template head<2>(), positions.template tail<2>());
    return g * deriv2 * g.transpose() + H * deriv1;
}

template <>
auto ESPCollisionTemplate<Vertex2, Edge2P1>::hessian(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params) const
    -> MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>
{
    assert(
        point_edge_distance_type_exact(
            positions.template head<2>(), positions.template segment<2>(2),
            positions.template segment<2>(4))
        == PointEdgeDistanceType::P_E);
    // Vertex2-Edge2P1 is constructed only at interior P_E.
    constexpr auto dtype = PointEdgeDistanceType::P_E;
    const double dist = std::sqrt(point_edge_distance(
        positions.template head<2>(), positions.template segment<2>(2),
        positions.template segment<2>(4), dtype));
    const double eps = params.dhat;
    double deriv1 = params.barrier->first_derivative(dist, eps);
    double deriv2 = params.barrier->second_derivative(dist, eps);
    deriv2 = deriv2 / (4 * dist * dist) - deriv1 / (4 * dist * dist * dist);
    deriv1 /= (2 * dist);
    const VectorMax9d g = point_edge_distance_gradient(
        positions.template head<2>(), positions.template segment<2>(2),
        positions.template segment<2>(4), dtype);
    const MatrixMax9d H = point_edge_distance_hessian(
        positions.template head<2>(), positions.template segment<2>(2),
        positions.template segment<2>(4), dtype);
    return g * deriv2 * g.transpose() + H * deriv1;
}

// ---- explicit instantiations ----

template class ESPCollisionTemplate<Vertex3, Vertex3>;
template class ESPCollisionTemplate<Edge3P1, Vertex3>;
template class ESPCollisionTemplate<Face3P1, Vertex3>;
template class ESPCollisionTemplate<Vertex2, Vertex2>;
template class ESPCollisionTemplate<Vertex2, Edge2P1>;

} // namespace ipc
