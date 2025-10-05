#include "smooth_collision.hpp"

namespace ipc {

template <typename PrimitiveA, typename PrimitiveB>
CollisionType SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::type() const
{
    if constexpr (
        std::is_same_v<PrimitiveA, Edge2> && std::is_same_v<PrimitiveB, Point2>)
        return CollisionType::EdgeVertex;
    if constexpr (
        std::is_same_v<PrimitiveA, Point2>
        && std::is_same_v<PrimitiveB, Point2>)
        return CollisionType::VertexVertex;
    if constexpr (
        std::is_same_v<PrimitiveA, Face> && std::is_same_v<PrimitiveB, Point3>)
        return CollisionType::FaceVertex;
    if constexpr (
        std::is_same_v<PrimitiveA, Edge3> && std::is_same_v<PrimitiveB, Point3>)
        return CollisionType::EdgeVertex;
    if constexpr (
        std::is_same_v<PrimitiveA, Edge3> && std::is_same_v<PrimitiveB, Edge3>)
        return CollisionType::EdgeEdge;
    if constexpr (
        std::is_same_v<PrimitiveA, Point3>
        && std::is_same_v<PrimitiveB, Point3>)
        return CollisionType::VertexVertex;

    throw std::runtime_error("Invalid collision pair type!");
    return CollisionType::VertexVertex;
}

Eigen::VectorXd SmoothCollision::dof(Eigen::ConstRef<Eigen::MatrixXd> X) const
{
    const int dim = X.cols();
    Eigen::VectorXd x(num_vertices() * dim);
    if (dim == 2) {
        for (int i = 0; i < num_vertices(); i++) {
            x.segment<2>(i * 2) = X.row(vertex_ids_[i]);
        }
    } else if (dim == 3) {
        for (int i = 0; i < num_vertices(); i++) {
            x.segment<3>(i * 3) = X.row(vertex_ids_[i]);
        }
    } else
        throw std::runtime_error("Invalid dimension!");
    return x;
}

template <typename PrimitiveA, typename PrimitiveB>
std::string SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::name() const
{
    if constexpr (
        std::is_same_v<PrimitiveA, Edge2> && std::is_same_v<PrimitiveB, Point2>) {
        return "edge-vert";
    }
    if constexpr (
        std::is_same_v<PrimitiveA, Point2>
        && std::is_same_v<PrimitiveB, Point2>) {
        return "vert-vert";
    }
    if constexpr (
        std::is_same_v<PrimitiveA, Face> && std::is_same_v<PrimitiveB, Point3>) {
        return "face-vert";
    }
    if constexpr (
        std::is_same_v<PrimitiveA, Edge3> && std::is_same_v<PrimitiveB, Point3>) {
        return "edge-vert";
    }
    if constexpr (
        std::is_same_v<PrimitiveA, Edge3> && std::is_same_v<PrimitiveB, Edge3>) {
        return "edge-edge";
    }
    if constexpr (
        std::is_same_v<PrimitiveA, Point3>
        && std::is_same_v<PrimitiveB, Point3>) {
        return "vert-vert";
    }

    throw std::runtime_error("Invalid collision pair type!");
    return "vert-vert";
}

template <typename PrimitiveA, typename PrimitiveB>
auto SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::get_core_indices() const
    -> Vector<int, n_core_dofs>
{
    Vector<int, n_core_dofs> core_indices;
    core_indices << Eigen::VectorXi::LinSpaced(
        n_core_dofs_A, 0, n_core_dofs_A - 1),
        Eigen::VectorXi::LinSpaced(
            n_core_dofs_B, pA->n_dofs(), pA->n_dofs() + n_core_dofs_B - 1);
    return core_indices;
}

template <typename PrimitiveA, typename PrimitiveB>
SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::SmoothCollisionTemplate(
    index_t primitive0_,
    index_t primitive1_,
    SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::DTYPE dtype,
    const CollisionMesh& mesh,
    const ParameterType& param,
    const double& dhat,
    const Eigen::MatrixXd& V)
    : SmoothCollision(primitive0_, primitive1_, dhat, mesh)
{
    VectorMax3d d =
        PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_closest_direction(
            mesh, V, primitive0_, primitive1_, dtype);
    pA = std::make_unique<PrimitiveA>(primitive0_, mesh, V, d, param);
    pB = std::make_unique<PrimitiveB>(primitive1_, mesh, V, -d, param);

    if ((pA->n_vertices() + pB->n_vertices()) * dim > ELEMENT_SIZE) {
        logger().error(
            "Too many neighbors for collision pair! {} > {}! Increase MAX_VERT_3D in common.hpp",
            pA->n_vertices() + pB->n_vertices(), MAX_VERT_3D);
    }

    int i = 0;
    Super::vertex_ids_.assign(
        pA->vertex_ids().size() + pB->vertex_ids().size(), -1);
    for (auto& v : pA->vertex_ids())
        Super::vertex_ids_[i++] = v;
    for (auto& v : pB->vertex_ids())
        Super::vertex_ids_[i++] = v;
    assert(i == pA->n_vertices() + pB->n_vertices());
    Super::is_active_ =
        (d.norm() < Super::dhat()) && pA->is_active() && pB->is_active();

    if (d.norm() < 1e-12) {
        logger().warn(
            "pair distance {}, id {} and {}, dtype {}, active {}", d.norm(),
            primitive0_, primitive1_,
            PrimitiveDistType<PrimitiveA, PrimitiveB>::name, Super::is_active_);

        logger().warn("value {}", (*this)(this->dof(V), param));
    }
}

template <typename PrimitiveA, typename PrimitiveB>
double SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::operator()(
    Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
    const ParameterType& params) const
{
    Vector<double, n_core_points * dim> x;
    x << positions.head(PrimitiveA::n_core_points * dim),
        positions.segment(pA->n_dofs(), PrimitiveB::n_core_points * dim);

    // grad of "d" wrt. points
    const Vector<double, dim> closest_direction =
        PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, double>::
            compute_closest_direction(x, DTYPE::AUTO);
    const double dist = closest_direction.norm();

    assert(positions.size() == pA->n_dofs() + pB->n_dofs());
    double a1 = pA->potential(closest_direction, positions.head(pA->n_dofs()));
    double a2 = pB->potential(-closest_direction, positions.tail(pB->n_dofs()));
    double a3 = Math<double>::inv_barrier(dist / Super::dhat(), params.r);
    double a4 =
        PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, double>::mollifier(
            x, dist * dist);

    if (params.r == 0) {
        logger().error("Invalid param!");
    }

    if (dist < 1e-12) {
        logger().warn(
            "pair distance {:.3e}, dhat {:.3e}, r {}, barrier {:.3e}, mollifier {:.3e}, orient {:.3e} {:.3e}",
            dist, Super::dhat(), params.r, a3, a4, a1, a2);
    }

    return a1 * a2 * a3 * a4;
}

template <typename PrimitiveA, typename PrimitiveB>
auto SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::gradient(
    Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
    const ParameterType& params) const -> Vector<double, -1, ELEMENT_SIZE>
{
    const auto core_indices = get_core_indices();

    Vector<double, n_core_dofs> x;
    x = positions(core_indices);

    const auto dtype =
        PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_distance_type(x);

    Vector<double, dim> closest_direction;
    Eigen::Matrix<double, dim, n_core_dofs> closest_direction_grad;
    std::tie(closest_direction, closest_direction_grad) = PrimitiveDistance<
        PrimitiveA, PrimitiveB>::compute_closest_direction_gradient(x, dtype);

    const double dist = closest_direction.norm();
    assert(dist > 0);

    // these two use autodiff with different variable count
    auto gA_reduced = pA->grad(closest_direction, positions.head(pA->n_dofs()));
    auto gB_reduced =
        pB->grad(-closest_direction, positions.tail(pB->n_dofs()));

    // gradient of barrier potential
    double barrier = 0;
    Vector<double, n_core_dofs> gBarrier = Vector<double, n_core_dofs>::Zero();
    {
        barrier = Math<double>::inv_barrier(dist / Super::dhat(), params.r);

        const Vector<double, dim> closest_direction_normalized =
            closest_direction / dist;
        const double barrier_1st_deriv =
            Math<double>::inv_barrier_grad(dist / Super::dhat(), params.r)
            / Super::dhat();
        const Vector<double, dim> gBarrier_wrt_d =
            barrier_1st_deriv * closest_direction_normalized;
        gBarrier = closest_direction_grad.transpose() * gBarrier_wrt_d;
    }

    // gradient of mollifier
    {
        double mollifier = 0;
        Vector<double, n_core_dofs> gMollifier =
            Vector<double, n_core_dofs>::Zero();
#ifdef DERIVATIVES_WITH_AUTODIFF
        DiffScalarBase::setVariableCount(n_core_dofs);
        using T = ADGrad<n_core_dofs>;
        Vector<T, n_core_dofs> xAD = slice_positions<T, n_core_dofs, 1>(x);
        Vector<T, dim> closest_direction_autodiff = PrimitiveDistanceTemplate<
            PrimitiveA, PrimitiveB, T>::compute_closest_direction(xAD, dtype);
        const auto dist_sqr_AD = closest_direction_autodiff.squaredNorm();
        auto mollifier_autodiff =
            PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::mollifier(
                xAD, dist_sqr_AD);
        mollifier = mollifier_autodiff.getValue();
        gMollifier = mollifier_autodiff.getGradient();
#else
        Vector<double, n_core_dofs + 1> mollifier_grad;
        std::tie(mollifier, mollifier_grad) = PrimitiveDistance<
            PrimitiveA, PrimitiveB>::compute_mollifier_gradient(x, dist * dist);

        const Vector<double, n_core_dofs> dist_sqr_grad =
            2 * closest_direction_grad.transpose() * closest_direction;
        mollifier_grad.head(n_core_dofs) +=
            mollifier_grad(n_core_dofs) * dist_sqr_grad;
        gMollifier = mollifier_grad.head(n_core_dofs);
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
                closest_direction_grad.transpose() * gA_reduced.head(dim);
            gA.head(pA->n_dofs()) += gA_reduced.tail(pA->n_dofs());

            gB(core_indices) =
                closest_direction_grad.transpose() * -gB_reduced.head(dim);
            gB.tail(pB->n_dofs()) += gB_reduced.tail(pB->n_dofs());
        }
        const double potential_a =
            pA->potential(closest_direction, positions.head(pA->n_dofs()));
        const double potential_b =
            pB->potential(-closest_direction, positions.tail(pB->n_dofs()));

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
    const ParameterType& params) const
    -> MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>
{
    const auto core_indices = get_core_indices();

    Vector<double, n_core_dofs> x;
    x = positions(core_indices);

    const auto dtype =
        PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_distance_type(x);

    Vector<double, dim> closest_direction;
    Eigen::Matrix<double, dim, n_core_dofs> closest_direction_grad;
    std::array<Eigen::Matrix<double, n_core_dofs, n_core_dofs>, dim>
        closest_direction_hess;
    std::tie(
        closest_direction, closest_direction_grad, closest_direction_hess) =
        PrimitiveDistance<PrimitiveA, PrimitiveB>::
            compute_closest_direction_hessian(x, dtype);

    const double dist = closest_direction.norm();

    // these two use autodiff with different variable count
    auto gA_reduced = pA->grad(closest_direction, positions.head(pA->n_dofs()));
    auto hA_reduced =
        pA->hessian(closest_direction, positions.head(pA->n_dofs()));
    auto gB_reduced =
        pB->grad(-closest_direction, positions.tail(pB->n_dofs()));
    auto hB_reduced =
        pB->hessian(-closest_direction, positions.tail(pB->n_dofs()));

    // hessian of barrier potential
    double barrier = 0;
    Vector<double, n_core_dofs> gBarrier = Vector<double, n_core_dofs>::Zero();
    Eigen::Matrix<double, n_core_dofs, n_core_dofs> hBarrier =
        Eigen::Matrix<double, n_core_dofs, n_core_dofs>::Zero();
    {
        barrier = Math<double>::inv_barrier(dist / Super::dhat(), params.r);

        const Vector<double, dim> closest_direction_normalized =
            closest_direction / dist;
        const double barrier_1st_deriv =
            Math<double>::inv_barrier_grad(dist / Super::dhat(), params.r)
            / Super::dhat();
        const Vector<double, dim> gBarrier_wrt_d =
            barrier_1st_deriv * closest_direction_normalized;
        gBarrier = closest_direction_grad.transpose() * gBarrier_wrt_d;

        const double barrier_2nd_deriv =
            Math<double>::inv_barrier_hess(dist / Super::dhat(), params.r)
            / Super::dhat() / Super::dhat();
        const Eigen::Matrix<double, dim, dim> hBarrier_wrt_d =
            (barrier_1st_deriv / dist)
                * Eigen::Matrix<double, dim, dim>::Identity()
            + (barrier_2nd_deriv - barrier_1st_deriv / dist)
                * closest_direction_normalized
                * closest_direction_normalized.transpose();
        hBarrier = closest_direction_grad.transpose() * hBarrier_wrt_d
            * closest_direction_grad;
        for (int d = 0; d < dim; d++)
            hBarrier += closest_direction_hess[d] * gBarrier_wrt_d(d);
    }

    // hessian of mollifier
    {
        double mollifier = 0;
        Vector<double, n_core_dofs> gMollifier =
            Vector<double, n_core_dofs>::Zero();
        Eigen::Matrix<double, n_core_dofs, n_core_dofs> hMollifier =
            Eigen::Matrix<double, n_core_dofs, n_core_dofs>::Zero();
#ifdef DERIVATIVES_WITH_AUTODIFF
        DiffScalarBase::setVariableCount(n_core_dofs);
        using T = ADHessian<n_core_dofs>;
        Vector<T, n_core_dofs> xAD = slice_positions<T, n_core_dofs, 1>(x);
        Vector<T, dim> closest_direction_autodiff = PrimitiveDistanceTemplate<
            PrimitiveA, PrimitiveB, T>::compute_closest_direction(xAD, dtype);
        const auto dist_sqr_AD = closest_direction_autodiff.squaredNorm();
        auto mollifier_autodiff =
            PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::mollifier(
                xAD, dist_sqr_AD);
        mollifier = mollifier_autodiff.getValue();

        gMollifier = mollifier_autodiff.getGradient();
        hMollifier = mollifier_autodiff.getHessian();
#else
        Vector<double, n_core_dofs + 1> mollifier_grad;
        Eigen::Matrix<double, n_core_dofs + 1, n_core_dofs + 1> mollifier_hess;
        std::tie(mollifier, mollifier_grad, mollifier_hess) = PrimitiveDistance<
            PrimitiveA, PrimitiveB>::compute_mollifier_hessian(x, dist * dist);

        const Vector<double, n_core_dofs> dist_sqr_grad =
            2 * closest_direction_grad.transpose() * closest_direction;
        mollifier_grad.head(n_core_dofs) +=
            mollifier_grad(n_core_dofs) * dist_sqr_grad;
        Eigen::Matrix<double, n_core_dofs, n_core_dofs> dist_sqr_hess =
            2 * closest_direction_grad.transpose() * closest_direction_grad;
        for (int d = 0; d < dim; d++)
            dist_sqr_hess +=
                2 * closest_direction(d) * closest_direction_hess[d];
        mollifier_hess.topLeftCorner(n_core_dofs, n_core_dofs) +=
            dist_sqr_hess * mollifier_grad(n_core_dofs)
            + dist_sqr_grad * mollifier_hess(n_core_dofs, n_core_dofs)
                * dist_sqr_grad.transpose()
            + dist_sqr_grad
                * mollifier_hess.block(n_core_dofs, 0, 1, n_core_dofs)
            + mollifier_hess.block(0, n_core_dofs, n_core_dofs, 1)
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
                closest_direction_grad.transpose() * gA_reduced.head(dim);
            gA.head(pA->n_dofs()) += gA_reduced.tail(pA->n_dofs());

            hA(core_indices, core_indices) = closest_direction_grad.transpose()
                * hA_reduced.topLeftCorner(dim, dim) * closest_direction_grad;
            for (int d = 0; d < dim; d++)
                hA(core_indices, core_indices) +=
                    gA_reduced(d) * closest_direction_hess[d];

            hA.topLeftCorner(pA->n_dofs(), pA->n_dofs()) +=
                hA_reduced.bottomRightCorner(pA->n_dofs(), pA->n_dofs());
            hA(core_indices, Eigen::seqN(0, pA->n_dofs())) +=
                closest_direction_grad.transpose()
                * hA_reduced.topRightCorner(dim, pA->n_dofs());
            hA(Eigen::seqN(0, pA->n_dofs()), core_indices) +=
                hA_reduced.bottomLeftCorner(pA->n_dofs(), dim)
                * closest_direction_grad;

            gB(core_indices) =
                closest_direction_grad.transpose() * -gB_reduced.head(dim);
            gB.tail(pB->n_dofs()) += gB_reduced.tail(pB->n_dofs());

            hB(core_indices, core_indices) = closest_direction_grad.transpose()
                * hB_reduced.topLeftCorner(dim, dim) * closest_direction_grad;
            for (int d = 0; d < dim; d++)
                hB(core_indices, core_indices) -=
                    gB_reduced(d) * closest_direction_hess[d];

            hB.bottomRightCorner(pB->n_dofs(), pB->n_dofs()) +=
                hB_reduced.bottomRightCorner(pB->n_dofs(), pB->n_dofs());
            hB(core_indices, Eigen::seqN(pA->n_dofs(), pB->n_dofs())) -=
                closest_direction_grad.transpose()
                * hB_reduced.topRightCorner(dim, pB->n_dofs());
            hB(Eigen::seqN(pA->n_dofs(), pB->n_dofs()), core_indices) -=
                hB_reduced.bottomLeftCorner(pB->n_dofs(), dim)
                * closest_direction_grad;
        }
        const double potential_a =
            pA->potential(closest_direction, positions.head(pA->n_dofs()));
        const double potential_b =
            pB->potential(-closest_direction, positions.tail(pB->n_dofs()));

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

    Vector<double, n_core_points * dim> x;
    x << positions.head(PrimitiveA::n_core_points * dim),
        positions.segment(pA->n_dofs(), PrimitiveB::n_core_points * dim);

    // grad of "d" wrt. points
    Vector<double, dim> closest_direction =
        PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, double>::
            compute_closest_direction(x, DTYPE::AUTO);

    return closest_direction.squaredNorm();
}

template <typename PrimitiveA, typename PrimitiveB>
auto SmoothCollisionTemplate<PrimitiveA, PrimitiveB>::core_vertex_ids() const
    -> std::array<index_t, n_core_dofs>
{
    std::array<index_t, n_core_dofs> vids;
    auto ids = get_core_indices();
    for (int i = 0; i < n_core_dofs; i++)
        vids[i] = Super::vertex_ids_[ids[i]];
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