#include "esp_collision_template.hpp"

#include <ipc/distance/distance_type_exact.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/esp/smooth_clamp.hpp>
#include <ipc/gcp/distance/point_edge.hpp>
#include <ipc/tangent/closest_point.hpp>
#include <ipc/utils/autodiff_types.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <algorithm>

namespace {

template <typename T> double scalar_val(const T& x)
{
    if constexpr (std::is_same_v<T, double>)
        return x;
    else
        return x.val;
}

// Evaluate barrier with AD or double types.
// NormalizedClampedLogBarrier must be checked before ClampedLogBarrier
// because the former inherits from the latter.
template <typename T>
T eval_barrier_ad(const ipc::Barrier& b, const T& dist, const T& dhat)
{
    using ipc::ClampedLogBarrier;
    using ipc::InversePowerBarrier;
    using ipc::NormalizedClampedLogBarrier;

    if (scalar_val(dist) >= scalar_val(dhat))
        return T(0.0);

    if (dynamic_cast<const NormalizedClampedLogBarrier<>*>(&b)) {
        const T t = dist / dhat;
        return -(t - 1.0) * (t - 1.0) * log(t);
    }
    if (dynamic_cast<const ClampedLogBarrier<>*>(&b)) {
        return -(dist - dhat) * (dist - dhat) * log(dist / dhat);
    }
    if (const auto* ipb = dynamic_cast<const InversePowerBarrier*>(&b)) {
        const double p = ipb->power();
        const T t = 2.0 * dist / dhat;
        T h;
        if (scalar_val(t) < 1.0) {
            h = 2.0 / 3.0 - t * t + t * t * t * 0.5;
        } else if (scalar_val(t) < 2.0) {
            const T s = 2.0 - t;
            h = s * s * s / 6.0;
        } else {
            return T(0.0);
        }
        return h * pow(dist, -p);
    }
    throw std::runtime_error("eval_barrier_ad: unsupported barrier type");
}

// Edge-Vertex 3D energy with AD types.
// positions order: [e0 (0:3), e1 (3:6), vertex (6:9)]
template <typename T>
T eval_ev3d_energy_ad(
    Eigen::ConstRef<ipc::VectorMax<double, ipc::ESPCollision::ELEMENT_SIZE>>
        positions,
    const ipc::ESPParameters& params,
    const ipc::AdaptiveSupport& adaptive,
    ipc::index_t edge_id)
{
    using Vec3T = Eigen::Vector3<T>;
    ipc::ScalarBase::setVariableCount(9);

    Vec3T e0, e1, p;
    for (int i = 0; i < 3; i++) {
        e0[i] = T(positions[i], i);
        e1[i] = T(positions[3 + i], 3 + i);
        p[i] = T(positions[6 + i], 6 + i);
    }

    // ESPCollisionTemplate<Edge3P1, Vertex3> is constructed only when
    // the closest point is in the interior of the edge (P_E in
    // ESPCollisionsBuilder<3>::reduce_point_edge_collision); endpoint
    // cases are reduced to Vertex3-Vertex3. So we always use the interior
    // projection here.
    const Vec3T t_edge = e1 - e0;
    const T u_raw = (p - e0).dot(t_edge) / t_edge.squaredNorm();
    const Vec3T closest = e0 + u_raw * t_edge;

    const T dist = sqrt((p - closest).squaredNorm());
    const T u_smooth = ipc::smooth_clamp01(u_raw);
    const T eps = (1.0 - u_smooth) * adaptive.edge(edge_id, 0.0)
        + u_smooth * adaptive.edge(edge_id, 1.0);

    params.record_dist(scalar_val(dist));
    return eval_barrier_ad(*params.barrier, dist, eps);
}

// Face-Vertex 3D energy with AD types.
// positions order: [f0 (0:3), f1 (3:6), f2 (6:9), vertex (9:12)]
template <typename T>
T eval_fv3d_energy_ad(
    Eigen::ConstRef<ipc::VectorMax<double, ipc::ESPCollision::ELEMENT_SIZE>>
        positions,
    const ipc::ESPParameters& params,
    const ipc::AdaptiveSupport& adaptive,
    ipc::index_t face_id)
{
    using Vec3T = Eigen::Vector3<T>;
    ipc::ScalarBase::setVariableCount(12);

    Vec3T f0, f1, f2, p;
    for (int i = 0; i < 3; i++) {
        f0[i] = T(positions[i], i);
        f1[i] = T(positions[3 + i], 3 + i);
        f2[i] = T(positions[6 + i], 6 + i);
        p[i] = T(positions[9 + i], 9 + i);
    }

    // ESPCollisionTemplate<Face3P1, Vertex3> is constructed only when
    // the closest point is in the interior of the triangle (P_T in
    // ESPCollisionsBuilder<3>::reduce_point_triangle_collision); edge
    // and vertex cases reduce to Edge3P1-Vertex3 / Vertex3-Vertex3. So we
    // always use the interior 2x2 solve here.
    const Vec3T e0t = f1 - f0, e1t = f2 - f0, dp = p - f0;
    const T A00 = e0t.dot(e0t), A01 = e0t.dot(e1t), A11 = e1t.dot(e1t);
    const T b0 = dp.dot(e0t), b1 = dp.dot(e1t);
    const T det = A00 * A11 - A01 * A01;
    const T u_raw = (b0 * A11 - b1 * A01) / det;
    const T v_raw = (b1 * A00 - b0 * A01) / det;
    const Vec3T closest = f0 + u_raw * e0t + v_raw * e1t;

    T u, v;
    ipc::smooth_clamp_simplex(u_raw, v_raw, u, v);

    const T dist = sqrt((p - closest).squaredNorm());
    const T eps = (1.0 - u - v) * adaptive.face(face_id, 0.0, 0.0)
        + u * adaptive.face(face_id, 1.0, 0.0)
        + v * adaptive.face(face_id, 0.0, 1.0);

    params.record_dist(scalar_val(dist));
    return eval_barrier_ad(*params.barrier, dist, eps);
}

// Vertex-Edge 2D energy with AD types.
// positions order: [q (0:2), e0 (2:4), e1 (4:6)]
template <typename T>
T eval_ve2d_energy_ad(
    Eigen::ConstRef<ipc::VectorMax<double, ipc::ESPCollision::ELEMENT_SIZE>>
        positions,
    const ipc::ESPParameters& params,
    const ipc::AdaptiveSupport& adaptive,
    ipc::index_t edge_id)
{
    using Vec2T = Eigen::Vector2<T>;
    ipc::ScalarBase::setVariableCount(6);

    Vec2T q, e0, e1;
    for (int i = 0; i < 2; i++) {
        q[i] = T(positions[i], i);
        e0[i] = T(positions[2 + i], 2 + i);
        e1[i] = T(positions[4 + i], 4 + i);
    }

    // ESPCollisionTemplate<Vertex2, Edge2P1> is constructed only when
    // the closest point is in the interior of the edge (the 2D edge-QP
    // builder in quadrature_potential.cpp routes endpoint cases to
    // Vertex2-Vertex2). So we always use the interior projection here.
    const Vec2T t_edge = e1 - e0;
    const T u_raw = (q - e0).dot(t_edge) / t_edge.squaredNorm();
    const Vec2T closest = e0 + u_raw * t_edge;

    const T dist = sqrt((q - closest).squaredNorm());
    const T u_smooth = ipc::smooth_clamp01(u_raw);
    const T eps = (1.0 - u_smooth) * adaptive.edge(edge_id, 0.0)
        + u_smooth * adaptive.edge(edge_id, 1.0);

    params.record_dist(scalar_val(dist));
    return eval_barrier_ad(*params.barrier, dist, eps);
}

} // anonymous namespace

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
    : primitive_a(_primitive0, mesh)
    , primitive_b(_primitive1, mesh)
{
    static_assert(
        Eigen::internal::packet_traits<double>::size == 1,
        "Eigen vectorization is NOT disabled!");
}

template <>
ESPCollisionTemplate<Vertex3, Vertex3>::ESPCollisionTemplate(
    index_t _primitive0, index_t _primitive1, const CollisionMesh& mesh)
    : primitive_a(std::min(_primitive0, _primitive1), mesh)
    , primitive_b(std::max(_primitive0, _primitive1), mesh)
{
}

// ---- vertex_id ----

template <typename PrimitiveA, typename PrimitiveB>
index_t ESPCollisionTemplate<PrimitiveA, PrimitiveB>::vertex_id(index_t i) const
{
    if (i < (index_t)primitive_a.n_vertices()) {
        return primitive_a.vertex_ids()[i];
    }
    i -= primitive_a.n_vertices();
    assert((index_t)primitive_b.n_vertices() > i);
    return primitive_b.vertex_ids()[i];
}

// ---- generic stubs ----

template <typename PrimitiveA, typename PrimitiveB>
double ESPCollisionTemplate<PrimitiveA, PrimitiveB>::operator()(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> /*positions*/,
    const ESPParameters& /*params*/,
    const AdaptiveSupport* /*adaptive*/) const
{
    return 0;
}

template <typename PrimitiveA, typename PrimitiveB>
auto ESPCollisionTemplate<PrimitiveA, PrimitiveB>::gradient(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> /*positions*/,
    const ESPParameters& /*params*/,
    const AdaptiveSupport* /*adaptive*/) const
    -> VectorMax<double, ELEMENT_SIZE>
{
    return VectorMax<double, ELEMENT_SIZE>::Zero(n_dofs());
}

template <typename PrimitiveA, typename PrimitiveB>
auto ESPCollisionTemplate<PrimitiveA, PrimitiveB>::hessian(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> /*positions*/,
    const ESPParameters& /*params*/,
    const AdaptiveSupport* /*adaptive*/) const
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
double ESPCollisionTemplate<Edge3P1, Edge3P1>::compute_distance(
    Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    const int n_verts = vertices.rows();
    if (n_verts > vertex_id(0) && n_verts > vertex_id(1)
        && n_verts > vertex_id(2) && n_verts > vertex_id(3)) {
        return edge_edge_distance(
            vertices.row(vertex_id(0)), vertices.row(vertex_id(1)),
            vertices.row(vertex_id(2)), vertices.row(vertex_id(3)));
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
    const ESPParameters& params,
    const AdaptiveSupport* adaptive) const
{
    const double dist =
        (positions.template head<3>() - positions.template segment<3>(3))
            .norm();
    const double eps =
        adaptive ? adaptive->vertex(primitive_a.id()) : params.dhat;
    params.record_dist(dist);
    return (*params.barrier)(dist, eps);
}

template <>
double ESPCollisionTemplate<Edge3P1, Vertex3>::operator()(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const AdaptiveSupport* adaptive) const
{
    assert(
        point_edge_distance_type_exact(
            positions.template segment<3>(6), positions.template head<3>(),
            positions.template segment<3>(3))
        == PointEdgeDistanceType::P_E);
    double eps;
    if (adaptive) {
        const double u = smooth_clamp01(point_edge_closest_point(
            positions.template segment<3>(6), positions.template head<3>(),
            positions.template segment<3>(3)));
        eps = adaptive->edge(primitive_a.id(), u);
    } else
        eps = params.dhat;
    // Edge3P1-Vertex3 is constructed only at interior P_E (see
    // ESPCollisionsBuilder<3>::reduce_point_edge_collision).
    const double dist = sqrt(point_edge_distance(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3)));
    params.record_dist(dist);
    return (*params.barrier)(dist, eps);
}

template <>
double ESPCollisionTemplate<Face3P1, Vertex3>::operator()(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const AdaptiveSupport* adaptive) const
{
    assert(
        point_triangle_distance_type_exact(
            positions.template segment<3>(9), positions.template head<3>(),
            positions.template segment<3>(3), positions.template segment<3>(6))
        == PointTriangleDistanceType::P_T);
    double eps;
    if (adaptive) {
        const Eigen::Vector2d uv_raw = point_triangle_closest_point(
            positions.template segment<3>(9), positions.template head<3>(),
            positions.template segment<3>(3), positions.template segment<3>(6));
        double u, v;
        smooth_clamp_simplex(uv_raw[0], uv_raw[1], u, v);
        eps = adaptive->face(primitive_a.id(), u, v);
    } else
        eps = params.dhat;
    // Face3P1-Vertex3 is constructed only at interior P_T (see
    // ESPCollisionsBuilder<3>::reduce_point_triangle_collision).
    const double dist = sqrt(point_triangle_distance(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6)));
    params.record_dist(dist);
    return (*params.barrier)(dist, eps);
}

template <>
auto ESPCollisionTemplate<Vertex3, Vertex3>::gradient(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const AdaptiveSupport* adaptive) const -> VectorMax<double, ELEMENT_SIZE>
{
    assert(positions.size() == 6);
    const double dist =
        (positions.template head<3>() - positions.template tail<3>()).norm();
    const double eps =
        adaptive ? adaptive->vertex(primitive_a.id()) : params.dhat;
    params.record_dist(dist);
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
    const ESPParameters& params,
    const AdaptiveSupport* adaptive) const -> VectorMax<double, ELEMENT_SIZE>
{
    assert(positions.size() == 9);
    assert(
        point_edge_distance_type_exact(
            positions.template segment<3>(6), positions.template head<3>(),
            positions.template segment<3>(3))
        == PointEdgeDistanceType::P_E);
    if (adaptive) {
        ScalarBase::setVariableCount(9);
        using T = ADGrad<9>;
        const T energy = eval_ev3d_energy_ad<T>(
            positions, params, *adaptive, primitive_a.id());
        return energy.grad;
    }
    // Edge3P1-Vertex3 is constructed only at interior P_E.
    constexpr auto dtype = PointEdgeDistanceType::P_E;
    const double dist = sqrt(point_edge_distance(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3), dtype));
    const double eps = params.dhat;
    params.record_dist(dist);
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
    const ESPParameters& params,
    const AdaptiveSupport* adaptive) const -> VectorMax<double, ELEMENT_SIZE>
{
    assert(positions.size() == 12);
    assert(
        point_triangle_distance_type_exact(
            positions.template segment<3>(9), positions.template head<3>(),
            positions.template segment<3>(3), positions.template segment<3>(6))
        == PointTriangleDistanceType::P_T);
    if (adaptive) {
        ScalarBase::setVariableCount(12);
        using T = ADGrad<12>;
        const T energy = eval_fv3d_energy_ad<T>(
            positions, params, *adaptive, primitive_a.id());
        return energy.grad;
    }
    // Face3P1-Vertex3 is constructed only at interior P_T.
    constexpr auto dtype = PointTriangleDistanceType::P_T;
    const double dist = sqrt(point_triangle_distance(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6),
        dtype));
    const double eps = params.dhat;
    params.record_dist(dist);
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
    const ESPParameters& params,
    const AdaptiveSupport* adaptive) const
    -> MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>
{
    assert(positions.size() == 6);
    const double dist =
        (positions.template head<3>() - positions.template tail<3>()).norm();
    const double eps =
        adaptive ? adaptive->vertex(primitive_a.id()) : params.dhat;
    params.record_dist(dist);
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
    const ESPParameters& params,
    const AdaptiveSupport* adaptive) const
    -> MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>
{
    assert(positions.size() == 9);
    assert(
        point_edge_distance_type_exact(
            positions.template segment<3>(6), positions.template head<3>(),
            positions.template segment<3>(3))
        == PointEdgeDistanceType::P_E);
    if (adaptive) {
        ScalarBase::setVariableCount(9);
        using T = ADHessian<9>;
        const T energy = eval_ev3d_energy_ad<T>(
            positions, params, *adaptive, primitive_a.id());
        return energy.Hess;
    }
    // Edge3P1-Vertex3 is constructed only at interior P_E.
    constexpr auto dtype = PointEdgeDistanceType::P_E;
    const double dist = sqrt(point_edge_distance(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3), dtype));
    const double eps = params.dhat;
    params.record_dist(dist);
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
    const ESPParameters& params,
    const AdaptiveSupport* adaptive) const
    -> MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>
{
    assert(positions.size() == 12);
    assert(
        point_triangle_distance_type_exact(
            positions.template segment<3>(9), positions.template head<3>(),
            positions.template segment<3>(3), positions.template segment<3>(6))
        == PointTriangleDistanceType::P_T);
    if (adaptive) {
        ScalarBase::setVariableCount(12);
        using T = ADHessian<12>;
        const T energy = eval_fv3d_energy_ad<T>(
            positions, params, *adaptive, primitive_a.id());
        return energy.Hess;
    }
    // Face3P1-Vertex3 is constructed only at interior P_T.
    constexpr auto dtype = PointTriangleDistanceType::P_T;
    const double dist = sqrt(point_triangle_distance(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6),
        dtype));
    const double eps = params.dhat;
    params.record_dist(dist);
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
    const AdaptiveSupport* adaptive,
    const NearFarBarrier* nf_barrier) const
{
    const double dist =
        (positions.template head<3>() - positions.template segment<3>(3))
            .norm();
    const double eps =
        adaptive ? adaptive->vertex(primitive_a.id()) : params.dhat;
    params.record_dist(dist);
    return { nf_barrier->near(dist, eps), nf_barrier->far(dist, eps) };
}

template <>
std::pair<double, double>
ESPCollisionTemplate<Edge3P1, Vertex3>::operator_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const AdaptiveSupport* adaptive,
    const NearFarBarrier* nf_barrier) const
{
    const double dist = sqrt(point_edge_distance(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3)));
    const double eps =
        adaptive ? adaptive->edge(primitive_a.id(), 0.5) : params.dhat;
    params.record_dist(dist);
    return { nf_barrier->near(dist, eps), nf_barrier->far(dist, eps) };
}

template <>
std::pair<double, double>
ESPCollisionTemplate<Face3P1, Vertex3>::operator_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const AdaptiveSupport* adaptive,
    const NearFarBarrier* nf_barrier) const
{
    const double dist = sqrt(point_triangle_distance(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6)));
    const double eps = adaptive
        ? adaptive->face(primitive_a.id(), 1.0 / 3.0, 1.0 / 3.0)
        : params.dhat;
    params.record_dist(dist);
    return { nf_barrier->near(dist, eps), nf_barrier->far(dist, eps) };
}

template <>
std::pair<
    VectorMax<double, ESPCollision::ELEMENT_SIZE>,
    VectorMax<double, ESPCollision::ELEMENT_SIZE>>
ESPCollisionTemplate<Vertex3, Vertex3>::gradient_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const AdaptiveSupport* adaptive,
    const NearFarBarrier* nf_barrier) const
{
    assert(positions.size() == 6);
    const double dist =
        (positions.template head<3>() - positions.template tail<3>()).norm();
    const double eps =
        adaptive ? adaptive->vertex(primitive_a.id()) : params.dhat;
    params.record_dist(dist);
    const double deriv_near =
        nf_barrier->first_derivative_near(dist, eps) / (dist * 2.);
    const double deriv_far =
        nf_barrier->first_derivative_far(dist, eps) / (dist * 2.);
    Vector6d g = point_point_distance_gradient(
        positions.template head<3>(), positions.template tail<3>());
    VectorMax<double, ESPCollision::ELEMENT_SIZE> g_near(6), g_far(6);
    g_near.head(6) = deriv_near * g;
    g_far.head(6) = deriv_far * g;
    return { g_near, g_far };
}

template <>
std::pair<
    VectorMax<double, ESPCollision::ELEMENT_SIZE>,
    VectorMax<double, ESPCollision::ELEMENT_SIZE>>
ESPCollisionTemplate<Edge3P1, Vertex3>::gradient_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const AdaptiveSupport* adaptive,
    const NearFarBarrier* nf_barrier) const
{
    assert(positions.size() == 9);
    auto dtype = point_edge_distance_type_exact(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3));
    const double dist = sqrt(point_edge_distance(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3), dtype));
    const double eps =
        adaptive ? adaptive->edge(primitive_a.id(), 0.5) : params.dhat;
    params.record_dist(dist);
    const double deriv_near =
        nf_barrier->first_derivative_near(dist, eps) / (dist * 2.);
    const double deriv_far =
        nf_barrier->first_derivative_far(dist, eps) / (dist * 2.);
    Vector9d g = point_edge_distance_gradient(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3), dtype);
    Vector9d g_near = deriv_near * g;
    Vector9d g_far = deriv_far * g;
    g_near = g_near({ 3, 4, 5, 6, 7, 8, 0, 1, 2 }).eval();
    g_far = g_far({ 3, 4, 5, 6, 7, 8, 0, 1, 2 }).eval();
    VectorMax<double, ESPCollision::ELEMENT_SIZE> result_near(9), result_far(9);
    result_near.head(9) = g_near;
    result_far.head(9) = g_far;
    return { result_near, result_far };
}

template <>
std::pair<
    VectorMax<double, ESPCollision::ELEMENT_SIZE>,
    VectorMax<double, ESPCollision::ELEMENT_SIZE>>
ESPCollisionTemplate<Face3P1, Vertex3>::gradient_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const AdaptiveSupport* adaptive,
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
    const double eps = adaptive
        ? adaptive->face(primitive_a.id(), 1.0 / 3.0, 1.0 / 3.0)
        : params.dhat;
    params.record_dist(dist);
    const double deriv_near =
        nf_barrier->first_derivative_near(dist, eps) / (dist * 2.);
    const double deriv_far =
        nf_barrier->first_derivative_far(dist, eps) / (dist * 2.);
    Vector12d g = point_triangle_distance_gradient(
        positions.template segment<3>(9), positions.template head<3>(),
        positions.template segment<3>(3), positions.template segment<3>(6),
        dtype);
    Vector12d g_near = deriv_near * g;
    Vector12d g_far = deriv_far * g;
    g_near = g_near({ 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2 }).eval();
    g_far = g_far({ 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2 }).eval();
    VectorMax<double, ESPCollision::ELEMENT_SIZE> result_near(12),
        result_far(12);
    result_near.head(12) = g_near;
    result_far.head(12) = g_far;
    return { result_near, result_far };
}

template <>
std::pair<
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>,
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>>
ESPCollisionTemplate<Vertex3, Vertex3>::hessian_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const AdaptiveSupport* adaptive,
    const NearFarBarrier* nf_barrier) const
{
    assert(positions.size() == 6);
    const double dist =
        (positions.template head<3>() - positions.template tail<3>()).norm();
    const double eps =
        adaptive ? adaptive->vertex(primitive_a.id()) : params.dhat;
    params.record_dist(dist);
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
    Matrix6d hess_near = g * deriv2_near * g.transpose() + h * deriv1_near;
    Matrix6d hess_far = g * deriv2_far * g.transpose() + h * deriv1_far;
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>
        result_near(6, 6), result_far(6, 6);
    result_near.block<6, 6>(0, 0) = hess_near;
    result_far.block<6, 6>(0, 0) = hess_far;
    return { result_near, result_far };
}

template <>
std::pair<
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>,
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>>
ESPCollisionTemplate<Edge3P1, Vertex3>::hessian_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const AdaptiveSupport* adaptive,
    const NearFarBarrier* nf_barrier) const
{
    assert(positions.size() == 9);
    auto dtype = point_edge_distance_type_exact(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3));
    const double dist = sqrt(point_edge_distance(
        positions.template segment<3>(6), positions.template head<3>(),
        positions.template segment<3>(3), dtype));
    const double eps =
        adaptive ? adaptive->edge(primitive_a.id(), 0.5) : params.dhat;
    params.record_dist(dist);
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
    Matrix9d hess_near = g * deriv2_near * g.transpose() + h * deriv1_near;
    Matrix9d hess_far = g * deriv2_far * g.transpose() + h * deriv1_far;
    std::vector<int> reorder { 3, 4, 5, 6, 7, 8, 0, 1, 2 };
    hess_near = hess_near(reorder, reorder).eval();
    hess_far = hess_far(reorder, reorder).eval();
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>
        result_near(9, 9), result_far(9, 9);
    result_near.block<9, 9>(0, 0) = hess_near;
    result_far.block<9, 9>(0, 0) = hess_far;
    return { result_near, result_far };
}

template <>
std::pair<
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>,
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>>
ESPCollisionTemplate<Face3P1, Vertex3>::hessian_nearfar(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const AdaptiveSupport* adaptive,
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
    const double eps = adaptive
        ? adaptive->face(primitive_a.id(), 1.0 / 3.0, 1.0 / 3.0)
        : params.dhat;
    params.record_dist(dist);
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
    Matrix12d hess_near = g * deriv2_near * g.transpose() + h * deriv1_near;
    Matrix12d hess_far = g * deriv2_far * g.transpose() + h * deriv1_far;
    std::vector<int> reorder { 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2 };
    hess_near = hess_near(reorder, reorder).eval();
    hess_far = hess_far(reorder, reorder).eval();
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>
        result_near(12, 12), result_far(12, 12);
    result_near.block<12, 12>(0, 0) = hess_near;
    result_far.block<12, 12>(0, 0) = hess_far;
    return { result_near, result_far };
}

// ---- 2D specializations ----
// positions layout VV: [q_x, q_y, v_x, v_y]
// positions layout VE: [q_x, q_y, e0_x, e0_y, e1_x, e1_y]

template <>
double ESPCollisionTemplate<Vertex2, Vertex2>::compute_distance(
    Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    const int n = vertices.rows();
    if (vertex_id(0) >= n || vertex_id(1) >= n)
        return std::numeric_limits<double>::max();
    return point_point_distance(
        vertices.row(vertex_id(0)), vertices.row(vertex_id(1)));
}

template <>
double ESPCollisionTemplate<Vertex2, Edge2P1>::compute_distance(
    Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    const int n = vertices.rows();
    if (vertex_id(0) >= n || vertex_id(1) >= n || vertex_id(2) >= n)
        return std::numeric_limits<double>::max();
    return point_edge_distance(
        vertices.row(vertex_id(0)), vertices.row(vertex_id(1)),
        vertices.row(vertex_id(2)));
}

template <>
double ESPCollisionTemplate<Vertex2, Vertex2>::operator()(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const AdaptiveSupport* adaptive) const
{
    const double dist =
        (positions.template head<2>() - positions.template tail<2>()).norm();
    const double eps = adaptive ? adaptive->vertex(primitive_b.id())
                                : // TODO Check primitive index
        params.dhat;
    params.record_dist(dist);
    return (*params.barrier)(dist, eps);
}

template <>
double ESPCollisionTemplate<Vertex2, Edge2P1>::operator()(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const AdaptiveSupport* adaptive) const
{
    assert(
        point_edge_distance_type_exact(
            positions.template head<2>(), positions.template segment<2>(2),
            positions.template segment<2>(4))
        == PointEdgeDistanceType::P_E);
    double eps;
    if (adaptive) {
        const double u = smooth_clamp01(point_edge_closest_point(
            positions.template head<2>(), positions.template segment<2>(2),
            positions.template segment<2>(4)));
        eps = adaptive->edge(primitive_b.id(), u);
    } else
        eps = params.dhat;
    // Vertex2-Edge2P1 is constructed only at interior P_E (the 2D edge-QP
    // builder routes endpoint cases to Vertex2-Vertex2).
    const double dist = std::sqrt(point_edge_distance(
        positions.template head<2>(), positions.template segment<2>(2),
        positions.template segment<2>(4)));
    params.record_dist(dist);
    return (*params.barrier)(dist, eps);
}

template <>
auto ESPCollisionTemplate<Vertex2, Vertex2>::gradient(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const AdaptiveSupport* adaptive) const -> VectorMax<double, ELEMENT_SIZE>
{
    const double dist =
        (positions.template head<2>() - positions.template tail<2>()).norm();
    const double eps =
        adaptive ? adaptive->vertex(primitive_b.id()) : params.dhat;
    params.record_dist(dist);
    const double deriv =
        params.barrier->first_derivative(dist, eps) / (dist * 2.0);
    const VectorMax6d g = point_point_distance_gradient(
        positions.template head<2>(), positions.template tail<2>());
    return deriv * g;
}

template <>
auto ESPCollisionTemplate<Vertex2, Edge2P1>::gradient(
    Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
    const ESPParameters& params,
    const AdaptiveSupport* adaptive) const -> VectorMax<double, ELEMENT_SIZE>
{
    assert(
        point_edge_distance_type_exact(
            positions.template head<2>(), positions.template segment<2>(2),
            positions.template segment<2>(4))
        == PointEdgeDistanceType::P_E);
    if (adaptive) {
        ScalarBase::setVariableCount(6);
        using T = ADGrad<6>;
        const T energy = eval_ve2d_energy_ad<T>(
            positions, params, *adaptive, primitive_b.id());
        return energy.grad;
    }
    // Vertex2-Edge2P1 is constructed only at interior P_E.
    constexpr auto dtype = PointEdgeDistanceType::P_E;
    const double dist = std::sqrt(point_edge_distance(
        positions.template head<2>(), positions.template segment<2>(2),
        positions.template segment<2>(4), dtype));
    const double eps = params.dhat;
    params.record_dist(dist);
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
    const ESPParameters& params,
    const AdaptiveSupport* adaptive) const
    -> MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>
{
    const double dist =
        (positions.template head<2>() - positions.template tail<2>()).norm();
    const double eps =
        adaptive ? adaptive->vertex(primitive_b.id()) : params.dhat;
    params.record_dist(dist);
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
    const ESPParameters& params,
    const AdaptiveSupport* adaptive) const
    -> MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>
{
    assert(
        point_edge_distance_type_exact(
            positions.template head<2>(), positions.template segment<2>(2),
            positions.template segment<2>(4))
        == PointEdgeDistanceType::P_E);
    if (adaptive) {
        ScalarBase::setVariableCount(6);
        using T = ADHessian<6>;
        const T energy = eval_ve2d_energy_ad<T>(
            positions, params, *adaptive, primitive_b.id());
        return energy.Hess;
    }
    // Vertex2-Edge2P1 is constructed only at interior P_E.
    constexpr auto dtype = PointEdgeDistanceType::P_E;
    const double dist = std::sqrt(point_edge_distance(
        positions.template head<2>(), positions.template segment<2>(2),
        positions.template segment<2>(4), dtype));
    const double eps = params.dhat;
    params.record_dist(dist);
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
