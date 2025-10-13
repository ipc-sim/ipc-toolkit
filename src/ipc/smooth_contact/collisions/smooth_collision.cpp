#include "smooth_collision.hpp"

namespace ipc {

// clang-format off
template <> CollisionType SmoothCollisionTemplate<Point2, Point2>::type() const { return CollisionType::VERTEX_VERTEX; }
template <> CollisionType SmoothCollisionTemplate<Point3, Point3>::type() const { return CollisionType::VERTEX_VERTEX; }
template <> CollisionType SmoothCollisionTemplate<Edge2, Point2>::type() const { return CollisionType::EDGE_VERTEX; }
template <> CollisionType SmoothCollisionTemplate<Edge3, Point3>::type() const { return CollisionType::EDGE_VERTEX; }
template <> CollisionType SmoothCollisionTemplate<Face, Point3>::type() const { return CollisionType::FACE_VERTEX; }
template <> CollisionType SmoothCollisionTemplate<Edge3, Edge3>::type() const { return CollisionType::EDGE_EDGE; }
// clang-format on

// clang-format off
template <> std::string SmoothCollisionTemplate<Point2, Point2>::name() const { return "vert-vert"; }
template <> std::string SmoothCollisionTemplate<Point3, Point3>::name() const { return "vert-vert"; }
template <> std::string SmoothCollisionTemplate<Edge2, Point2>::name() const { return "edge-vert"; }
template <> std::string SmoothCollisionTemplate<Edge3, Point3>::name() const { return "edge-vert"; }
template <> std::string SmoothCollisionTemplate<Face, Point3>::name() const { return "face-vert"; }
template <> std::string SmoothCollisionTemplate<Edge3, Edge3>::name() const { return "edge-edge"; }
// clang-format on

Eigen::VectorXd SmoothCollision::dof(Eigen::ConstRef<Eigen::MatrixXd> X) const
{
    const int DIM = X.cols();
    Eigen::VectorXd x(num_vertices() * DIM);
    if (DIM == 2) {
        for (int i = 0; i < num_vertices(); i++) {
            x.segment<2>(i * 2) = X.row(m_vertex_ids[i]);
        }
    } else if (DIM == 3) {
        for (int i = 0; i < num_vertices(); i++) {
            x.segment<3>(i * 3) = X.row(m_vertex_ids[i]);
        }
    } else {
        throw std::runtime_error("Invalid dimension!");
    }
    return x;
}

template <typename PrimitiveA, typename PrimitiveB>
auto SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::get_core_indices() const
    -> Vector<int, N_CORE_DOFS>
{
    Vector<int, N_CORE_DOFS> core_indices;
    core_indices << Eigen::VectorXi::LinSpaced(
        N_CORE_DOFS_A, 0, N_CORE_DOFS_A - 1),
        Eigen::VectorXi::LinSpaced(
            N_CORE_DOFS_B, primitive_a->n_dofs(),
            primitive_a->n_dofs() + N_CORE_DOFS_B - 1);
    return core_indices;
}

template <typename PrimitiveA, typename PrimitiveB>
SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::SmoothCollisionTemplate(
    index_t _primitive0,
    index_t _primitive1,
    SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::DTYPE dtype,
    const CollisionMesh& mesh,
    const SmoothContactParameters& params,
    const double _dhat,
    const Eigen::MatrixXd& V)
    : SmoothCollision(_primitive0, _primitive1, _dhat, mesh)
{
    VectorMax3d d =
        PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_closest_direction(
            mesh, V, _primitive0, _primitive1, dtype);
    primitive_a = std::make_unique<PrimitiveA>(_primitive0, mesh, V, d, params);
    primitive_b =
        std::make_unique<PrimitiveB>(_primitive1, mesh, V, -d, params);

    if ((primitive_a->n_vertices() + primitive_b->n_vertices()) * DIM
        > ELEMENT_SIZE) {
        logger().error(
            "Too many neighbors for collision pair! {} > {}! Increase MAX_VERT_3D in common.hpp",
            primitive_a->n_vertices() + primitive_b->n_vertices(), MAX_VERT_3D);
    }

    int i = 0;
    m_vertex_ids.assign(
        primitive_a->vertex_ids().size() + primitive_b->vertex_ids().size(),
        -1);
    for (auto& v : primitive_a->vertex_ids()) {
        m_vertex_ids[i++] = v;
    }
    for (auto& v : primitive_b->vertex_ids()) {
        m_vertex_ids[i++] = v;
    }
    assert(i == primitive_a->n_vertices() + primitive_b->n_vertices());
    m_is_active = (d.norm() < m_dhat) && primitive_a->is_active()
        && primitive_b->is_active();

    if (d.norm() < 1e-12) {
        logger().warn(
            "pair distance {}, id {} and {}, dtype {}, active {}", d.norm(),
            _primitive0, _primitive1,
            PrimitiveDistType<PrimitiveA, PrimitiveB>::NAME, m_is_active);

        logger().warn("value {}", (*this)(this->dof(V), params));
    }
}

template <typename PrimitiveA, typename PrimitiveB>
double SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::operator()(
    Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
    const SmoothContactParameters& params) const
{
    Vector<double, N_CORE_POINTS * DIM> x;
    x << positions.head(PrimitiveA::N_CORE_POINTS * DIM),
        positions.segment(
            primitive_a->n_dofs(), PrimitiveB::N_CORE_POINTS * DIM);

    // grad of "d" wrt. points
    const Vector<double, DIM> closest_direction =
        PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, double>::
            compute_closest_direction(x, DTYPE::AUTO);
    const double dist = closest_direction.norm();

    assert(positions.size() == primitive_a->n_dofs() + primitive_b->n_dofs());
    double a1 = primitive_a->potential(
        closest_direction, positions.head(primitive_a->n_dofs()));
    double a2 = primitive_b->potential(
        -closest_direction, positions.tail(primitive_b->n_dofs()));
    double a3 = Math<double>::inv_barrier(dist / dhat(), params.r);
    double a4 =
        PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, double>::mollifier(
            x, dist * dist);

    if (params.r == 0) {
        logger().error("Invalid params!");
    }

    if (dist < 1e-12) {
        logger().warn(
            "pair distance {:.3e}, dhat {:.3e}, r {}, barrier {:.3e}, mollifier {:.3e}, orient {:.3e} {:.3e}",
            dist, dhat(), params.r, a3, a4, a1, a2);
    }

    return a1 * a2 * a3 * a4;
}

template <typename PrimitiveA, typename PrimitiveB>
auto SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::gradient(
    Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
    const SmoothContactParameters& params) const
    -> Vector<double, -1, ELEMENT_SIZE>
{
    const auto core_indices = get_core_indices();

    Vector<double, N_CORE_DOFS> x;
    x = positions(core_indices);

    const auto dtype =
        PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_distance_type(x);

    Vector<double, DIM> closest_direction;
    Eigen::Matrix<double, DIM, N_CORE_DOFS> closest_direction_grad;
    std::tie(closest_direction, closest_direction_grad) = PrimitiveDistance<
        PrimitiveA, PrimitiveB>::compute_closest_direction_gradient(x, dtype);

    const double dist = closest_direction.norm();
    assert(dist > 0);

    // these two use autodiff with different variable count
    auto gA_reduced = primitive_a->grad(
        closest_direction, positions.head(primitive_a->n_dofs()));
    auto gB_reduced = primitive_b->grad(
        -closest_direction, positions.tail(primitive_b->n_dofs()));

    // gradient of barrier potential
    double barrier = 0;
    Vector<double, N_CORE_DOFS> gBarrier = Vector<double, N_CORE_DOFS>::Zero();
    {
        barrier = Math<double>::inv_barrier(dist / dhat(), params.r);

        const Vector<double, DIM> closest_direction_normalized =
            closest_direction / dist;
        const double barrier_1st_deriv =
            Math<double>::inv_barrier_grad(dist / dhat(), params.r) / dhat();
        const Vector<double, DIM> gBarrier_wrt_d =
            barrier_1st_deriv * closest_direction_normalized;
        gBarrier = closest_direction_grad.transpose() * gBarrier_wrt_d;
    }

    // gradient of mollifier
    {
        double mollifier = 0;
        Vector<double, N_CORE_DOFS> gMollifier =
            Vector<double, N_CORE_DOFS>::Zero();
#ifdef DERIVATIVES_WITH_AUTODIFF
        ScalarBase::setVariableCount(N_CORE_DOFS);
        using T = ADGrad<N_CORE_DOFS>;
        Vector<T, N_CORE_DOFS> xAD = slice_positions<T, N_CORE_DOFS, 1>(x);
        Vector<T, DIM> closest_direction_autodiff = PrimitiveDistanceTemplate<
            PrimitiveA, PrimitiveB, T>::compute_closest_direction(xAD, dtype);
        const auto dist_sqr_AD = closest_direction_autodiff.squaredNorm();
        auto mollifier_autodiff =
            PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::mollifier(
                xAD, dist_sqr_AD);
        mollifier = mollifier_autodiff.val;
        gMollifier = mollifier_autodiff.grad;
#else
        Vector<double, N_CORE_DOFS + 1> mollifier_grad;
        std::tie(mollifier, mollifier_grad) = PrimitiveDistance<
            PrimitiveA, PrimitiveB>::compute_mollifier_gradient(x, dist * dist);

        const Vector<double, N_CORE_DOFS> dist_sqr_grad =
            2 * closest_direction_grad.transpose() * closest_direction;
        mollifier_grad.head(N_CORE_DOFS) +=
            mollifier_grad(N_CORE_DOFS) * dist_sqr_grad;
        gMollifier = mollifier_grad.head(N_CORE_DOFS);
#endif
        // merge mollifier into barrier
        gBarrier = gBarrier * mollifier + gMollifier * barrier;
        barrier *= mollifier;
    }

    // grad of tangent/normal terms
    double orient = 0;
    Vector<double, -1, ELEMENT_SIZE> gOrient;
    {
        Vector<double, -1, ELEMENT_SIZE>
            gA = Vector<double, -1, ELEMENT_SIZE>::Zero(n_dofs()),
            gB = Vector<double, -1, ELEMENT_SIZE>::Zero(n_dofs());
        {
            gA(core_indices) =
                closest_direction_grad.transpose() * gA_reduced.head(DIM);
            gA.head(primitive_a->n_dofs()) +=
                gA_reduced.tail(primitive_a->n_dofs());

            gB(core_indices) =
                closest_direction_grad.transpose() * -gB_reduced.head(DIM);
            gB.tail(primitive_b->n_dofs()) +=
                gB_reduced.tail(primitive_b->n_dofs());
        }
        const double potential_a = primitive_a->potential(
            closest_direction, positions.head(primitive_a->n_dofs()));
        const double potential_b = primitive_b->potential(
            -closest_direction, positions.tail(primitive_b->n_dofs()));

        orient = potential_a * potential_b;
        gOrient = gA * potential_b + gB * potential_a;
    }

    // merge barrier into orient
    gOrient *= barrier;
    gOrient(core_indices) += gBarrier * orient;
    orient *= barrier;

    return gOrient;
}

template <typename PrimitiveA, typename PrimitiveB>
auto SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::hessian(
    Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
    const SmoothContactParameters& params) const
    -> MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>
{
    const auto core_indices = get_core_indices();

    Vector<double, N_CORE_DOFS> x;
    x = positions(core_indices);

    const auto dtype =
        PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_distance_type(x);

    Vector<double, DIM> closest_direction;
    Eigen::Matrix<double, DIM, N_CORE_DOFS> closest_direction_grad;
    std::array<Eigen::Matrix<double, N_CORE_DOFS, N_CORE_DOFS>, DIM>
        closest_direction_hess;
    std::tie(
        closest_direction, closest_direction_grad, closest_direction_hess) =
        PrimitiveDistance<PrimitiveA, PrimitiveB>::
            compute_closest_direction_hessian(x, dtype);

    const double dist = closest_direction.norm();

    // these two use autodiff with different variable count
    auto gA_reduced = primitive_a->grad(
        closest_direction, positions.head(primitive_a->n_dofs()));
    auto hA_reduced = primitive_a->hessian(
        closest_direction, positions.head(primitive_a->n_dofs()));
    auto gB_reduced = primitive_b->grad(
        -closest_direction, positions.tail(primitive_b->n_dofs()));
    auto hB_reduced = primitive_b->hessian(
        -closest_direction, positions.tail(primitive_b->n_dofs()));

    // hessian of barrier potential
    double barrier = 0;
    Vector<double, N_CORE_DOFS> gBarrier = Vector<double, N_CORE_DOFS>::Zero();
    Eigen::Matrix<double, N_CORE_DOFS, N_CORE_DOFS> hBarrier =
        Eigen::Matrix<double, N_CORE_DOFS, N_CORE_DOFS>::Zero();
    {
        barrier = Math<double>::inv_barrier(dist / dhat(), params.r);

        const Vector<double, DIM> closest_direction_normalized =
            closest_direction / dist;
        const double barrier_1st_deriv =
            Math<double>::inv_barrier_grad(dist / dhat(), params.r) / dhat();
        const Vector<double, DIM> gBarrier_wrt_d =
            barrier_1st_deriv * closest_direction_normalized;
        gBarrier = closest_direction_grad.transpose() * gBarrier_wrt_d;

        const double barrier_2nd_deriv =
            Math<double>::inv_barrier_hess(dist / dhat(), params.r) / dhat()
            / dhat();
        const Eigen::Matrix<double, DIM, DIM> hBarrier_wrt_d =
            (barrier_1st_deriv / dist)
                * Eigen::Matrix<double, DIM, DIM>::Identity()
            + (barrier_2nd_deriv - barrier_1st_deriv / dist)
                * closest_direction_normalized
                * closest_direction_normalized.transpose();
        hBarrier = closest_direction_grad.transpose() * hBarrier_wrt_d
            * closest_direction_grad;
        for (int d = 0; d < DIM; d++)
            hBarrier += closest_direction_hess[d] * gBarrier_wrt_d(d);
    }

    // hessian of mollifier
    {
        double mollifier = 0;
        Vector<double, N_CORE_DOFS> gMollifier =
            Vector<double, N_CORE_DOFS>::Zero();
        Eigen::Matrix<double, N_CORE_DOFS, N_CORE_DOFS> hMollifier =
            Eigen::Matrix<double, N_CORE_DOFS, N_CORE_DOFS>::Zero();
#ifdef DERIVATIVES_WITH_AUTODIFF
        ScalarBase::setVariableCount(N_CORE_DOFS);
        using T = ADHessian<N_CORE_DOFS>;
        Vector<T, N_CORE_DOFS> xAD = slice_positions<T, N_CORE_DOFS, 1>(x);
        Vector<T, DIM> closest_direction_autodiff = PrimitiveDistanceTemplate<
            PrimitiveA, PrimitiveB, T>::compute_closest_direction(xAD, dtype);
        const auto dist_sqr_AD = closest_direction_autodiff.squaredNorm();
        auto mollifier_autodiff =
            PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::mollifier(
                xAD, dist_sqr_AD);
        mollifier = mollifier_autodiff.val;

        gMollifier = mollifier_autodiff.grad;
        hMollifier = mollifier_autodiff.Hess;
#else
        Vector<double, N_CORE_DOFS + 1> mollifier_grad;
        Eigen::Matrix<double, N_CORE_DOFS + 1, N_CORE_DOFS + 1> mollifier_hess;
        std::tie(mollifier, mollifier_grad, mollifier_hess) = PrimitiveDistance<
            PrimitiveA, PrimitiveB>::compute_mollifier_hessian(x, dist * dist);

        const Vector<double, N_CORE_DOFS> dist_sqr_grad =
            2 * closest_direction_grad.transpose() * closest_direction;
        mollifier_grad.head(N_CORE_DOFS) +=
            mollifier_grad(N_CORE_DOFS) * dist_sqr_grad;
        Eigen::Matrix<double, N_CORE_DOFS, N_CORE_DOFS> dist_sqr_hess =
            2 * closest_direction_grad.transpose() * closest_direction_grad;
        for (int d = 0; d < DIM; d++)
            dist_sqr_hess +=
                2 * closest_direction(d) * closest_direction_hess[d];
        mollifier_hess.topLeftCorner(N_CORE_DOFS, N_CORE_DOFS) +=
            dist_sqr_hess * mollifier_grad(N_CORE_DOFS)
            + dist_sqr_grad * mollifier_hess(N_CORE_DOFS, N_CORE_DOFS)
                * dist_sqr_grad.transpose()
            + dist_sqr_grad
                * mollifier_hess.block(N_CORE_DOFS, 0, 1, N_CORE_DOFS)
            + mollifier_hess.block(0, N_CORE_DOFS, N_CORE_DOFS, 1)
                * dist_sqr_grad.transpose();

        gMollifier = mollifier_grad.head(core_indices.size());
        hMollifier = mollifier_hess.topLeftCorner(
            core_indices.size(), core_indices.size());
#endif
        // merge mollifier into barrier
        hBarrier = hBarrier * mollifier + gBarrier * gMollifier.transpose()
            + gMollifier * gBarrier.transpose() + hMollifier * barrier;
        gBarrier = gBarrier * mollifier + gMollifier * barrier;
        barrier *= mollifier;
    }

    // grad of tangent/normal terms
    double orient = 0;
    Vector<double, -1, ELEMENT_SIZE> gOrient;
    MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE> hOrient;
    {
        Vector<double, -1, ELEMENT_SIZE>
            gA = Vector<double, -1, ELEMENT_SIZE>::Zero(n_dofs()),
            gB = Vector<double, -1, ELEMENT_SIZE>::Zero(n_dofs());
        MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>
            hA = MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>::Zero(
                n_dofs(), n_dofs()),
            hB = MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>::Zero(
                n_dofs(), n_dofs());
        {
            gA(core_indices) =
                closest_direction_grad.transpose() * gA_reduced.head(DIM);
            gA.head(primitive_a->n_dofs()) +=
                gA_reduced.tail(primitive_a->n_dofs());

            hA(core_indices, core_indices) = closest_direction_grad.transpose()
                * hA_reduced.topLeftCorner(DIM, DIM) * closest_direction_grad;
            for (int d = 0; d < DIM; d++)
                hA(core_indices, core_indices) +=
                    gA_reduced(d) * closest_direction_hess[d];

            hA.topLeftCorner(primitive_a->n_dofs(), primitive_a->n_dofs()) +=
                hA_reduced.bottomRightCorner(
                    primitive_a->n_dofs(), primitive_a->n_dofs());
            hA(core_indices, Eigen::seqN(0, primitive_a->n_dofs())) +=
                closest_direction_grad.transpose()
                * hA_reduced.topRightCorner(DIM, primitive_a->n_dofs());
            hA(Eigen::seqN(0, primitive_a->n_dofs()), core_indices) +=
                hA_reduced.bottomLeftCorner(primitive_a->n_dofs(), DIM)
                * closest_direction_grad;

            gB(core_indices) =
                closest_direction_grad.transpose() * -gB_reduced.head(DIM);
            gB.tail(primitive_b->n_dofs()) +=
                gB_reduced.tail(primitive_b->n_dofs());

            hB(core_indices, core_indices) = closest_direction_grad.transpose()
                * hB_reduced.topLeftCorner(DIM, DIM) * closest_direction_grad;
            for (int d = 0; d < DIM; d++)
                hB(core_indices, core_indices) -=
                    gB_reduced(d) * closest_direction_hess[d];

            hB.bottomRightCorner(
                primitive_b->n_dofs(), primitive_b->n_dofs()) +=
                hB_reduced.bottomRightCorner(
                    primitive_b->n_dofs(), primitive_b->n_dofs());
            hB(core_indices,
               Eigen::seqN(primitive_a->n_dofs(), primitive_b->n_dofs())) -=
                closest_direction_grad.transpose()
                * hB_reduced.topRightCorner(DIM, primitive_b->n_dofs());
            hB(Eigen::seqN(primitive_a->n_dofs(), primitive_b->n_dofs()),
               core_indices) -=
                hB_reduced.bottomLeftCorner(primitive_b->n_dofs(), DIM)
                * closest_direction_grad;
        }
        const double potential_a = primitive_a->potential(
            closest_direction, positions.head(primitive_a->n_dofs()));
        const double potential_b = primitive_b->potential(
            -closest_direction, positions.tail(primitive_b->n_dofs()));

        orient = potential_a * potential_b;
        gOrient = gA * potential_b + gB * potential_a;
        hOrient = (potential_a * hB + potential_b * hA)
            + (gA * gB.transpose() + gB * gA.transpose());
    }

    // merge barrier into orient
    hOrient *= barrier;
    hOrient(core_indices, core_indices) += hBarrier * orient;
    hOrient(Eigen::all, core_indices) += gOrient * gBarrier.transpose();
    hOrient(core_indices, Eigen::all) += gBarrier * gOrient.transpose();
    gOrient *= barrier;
    gOrient(core_indices) += gBarrier * orient;
    orient *= barrier;

    // return hOrient
    return (hOrient + hOrient.transpose()) / 2.;
}

// ---- distance ----

template <typename PrimitiveA, typename PrimitiveB>
double SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::compute_distance(
    Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    Vector<double, -1, ELEMENT_SIZE> positions = dof(vertices);

    Vector<double, N_CORE_POINTS * DIM> x;
    x << positions.head(PrimitiveA::N_CORE_POINTS * DIM),
        positions.segment(
            primitive_a->n_dofs(), PrimitiveB::N_CORE_POINTS * DIM);

    // grad of "d" wrt. points
    Vector<double, DIM> closest_direction =
        PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, double>::
            compute_closest_direction(x, DTYPE::AUTO);

    return closest_direction.squaredNorm();
}

template <typename PrimitiveA, typename PrimitiveB>
auto SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::core_vertex_ids() const
    -> std::array<index_t, N_CORE_DOFS>
{
    std::array<index_t, N_CORE_DOFS> vids;
    auto ids = get_core_indices();
    for (int i = 0; i < N_CORE_DOFS; i++) {
        vids[i] = m_vertex_ids[ids[i]];
    }
    return vids;
}

// Note: Primitive pair order cannot change
template class SmoothCollisionTemplate<Edge2, Point2>;
template class SmoothCollisionTemplate<Point2, Point2>;

template class SmoothCollisionTemplate<Edge3, Point3>;
template class SmoothCollisionTemplate<Edge3, Edge3>;
template class SmoothCollisionTemplate<Point3, Point3>;
template class SmoothCollisionTemplate<Face, Point3>;
} // namespace ipc