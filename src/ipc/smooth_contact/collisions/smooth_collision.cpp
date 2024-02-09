#include "smooth_collision.hpp"
#include <ipc/smooth_contact/primitives/point3.hpp>
#include <ipc/smooth_contact/primitives/edge.hpp>
#include <ipc/smooth_contact/primitives/face.hpp>

namespace ipc {
namespace {
    template <int max_vert, typename T, int dim, int max_dim = dim>
    Vector<T, dim, max_dim> reorder_vector(
        const SmoothCollision<max_vert>& a,
        const SmoothCollision<max_vert>& b,
        const Vector<T, dim, max_dim>& vec_a)
    {
        auto ids_a = a.vertices;
        auto ids_b = b.vertices;

        auto vec_b = vec_a;
        for (int i = 0; i < ids_a.size(); i++) {
            if (ids_a[i] < 0)
                break;
            for (int j = 0; j < ids_b.size(); j++) {
                if (ids_b[j] == ids_a[i]) {
                    vec_b.segment(j * 3, 3) = vec_a.segment(i * 3, 3);
                    break;
                }
                if (ids_b[j] < 0)
                    break;
            }
        }

        return vec_b;
    }

    template <int max_vert, typename T, int dim>
    MatrixMax<T, dim, dim> reorder_matrix(
        const SmoothCollision<max_vert>& a,
        const SmoothCollision<max_vert>& b,
        const MatrixMax<T, dim, dim>& mat_a)
    {
        auto ids_a = a.vertices;
        auto ids_b = b.vertices;

        std::map<int, int> id_map;
        for (int i = 0; i < ids_a.size(); i++) {
            if (ids_a[i] < 0)
                break;
            for (int j = 0; j < ids_b.size(); j++) {
                if (ids_b[j] == ids_a[i]) {
                    id_map[i] = j;
                    break;
                }
                if (ids_b[j] < 0)
                    break;
            }
        }

        auto mat_b = mat_a;
        for (int i = 0; i < ids_a.size(); i++)
            for (int j = 0; j < ids_a.size(); j++) {
                if (ids_a[i] < 0 || ids_a[j] < 0)
                    break;
                mat_b.block(id_map[i] * 3, id_map[j] * 3, 3, 3) =
                    mat_a.block(i * 3, j * 3, 3, 3);
            }

        return mat_b;
    }
} // namespace

template <int max_vert, typename PrimitiveA, typename PrimitiveB>
SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::
    SmoothCollisionTemplate(
        long primitive0_,
        long primitive1_,
        SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::DTYPE dtype,
        const CollisionMesh& mesh,
        const ParameterType& param,
        const double& dhat,
        const Eigen::MatrixXd& V)
    : SmoothCollision<max_vert>(primitive0_, primitive1_, dhat, mesh)
{
    VectorMax3d d =
        PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_closest_direction(
            mesh, V, primitive0_, primitive1_, dtype);
    pA = std::make_shared<PrimitiveA>(
        primitive0_, mesh, V, d, param.alpha, param.beta);
    pB = std::make_shared<PrimitiveB>(
        primitive1_, mesh, V, -d, param.alpha, param.beta);
    int i = 0;
    for (auto& v : pA->vertex_ids())
        Super::vertices[i++] = v;
    for (auto& v : pB->vertex_ids())
        Super::vertices[i++] = v;
    assert(i == pA->n_vertices() + pB->n_vertices());
    Super::is_active_ =
        (d.norm() < Super::get_dhat()) && pA->is_active() && pB->is_active();

    if (d.norm() < 1e-13)
        logger().warn("pair distance {}, id {} and {}, dtype {}, active {}", d.norm(), primitive0_,
                        primitive1_, PrimitiveDistType<PrimitiveA, PrimitiveB>::name, Super::is_active_);

    if (!Super::is_active())
        return;

    // test derivatives
    // std::shared_ptr<SmoothCollision<max_vert>> debug_ptr;
    // if constexpr (std::is_same<PrimitiveA, Edge3>::value &&
    // std::is_same<PrimitiveB, Edge3>::value)
    //     debug_ptr =
    //     std::make_shared<SmoothEdgeEdge3Collision>(Super::primitive0,
    //     Super::primitive1, mesh, param, dhat, V);
    // else if constexpr (std::is_same<PrimitiveA, Face>::value &&
    // std::is_same<PrimitiveB, Point3>::value)
    //     debug_ptr =
    //     std::make_shared<SmoothFaceVertexCollision>(Super::primitive0,
    //     Super::primitive1, mesh, param, dhat, V);
    // else if constexpr (std::is_same<PrimitiveA, Edge3>::value &&
    // std::is_same<PrimitiveB, Point3>::value)
    //     debug_ptr =
    //     std::make_shared<SmoothEdgeVertex3Collision>(Super::primitive0,
    //     Super::primitive1, mesh, param, dhat, V);
    // else if constexpr (std::is_same<PrimitiveA, Point3>::value &&
    // std::is_same<PrimitiveB, Point3>::value)
    //     debug_ptr =
    //     std::make_shared<SmoothVertexVertex3Collision>(Super::primitive0,
    //     Super::primitive1, mesh, param, dhat, V);

    // if (!debug_ptr)
    //     return;

    // double val = (*this)(this->dof(V, mesh.edges(), mesh.faces()), param);
    // double cc_val = (*debug_ptr)(debug_ptr->dof(V, mesh.edges(),
    // mesh.faces()), param);

    // if (std::abs(cc_val - val) > 1e-6 * std::min(std::abs(cc_val),
    // std::abs(val)))
    //     logger().error("Inconsistent potential: {} vs {}", cc_val, val);
}

template <int max_vert, typename PrimitiveA, typename PrimitiveB>
SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::
    ~SmoothCollisionTemplate()
{
}

template <int max_vert, typename PrimitiveA, typename PrimitiveB>
double SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::operator()(
    const Vector<
        double,
        -1,
        SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>&
        positions,
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
    double a3 = Math<double>::inv_barrier(dist / Super::get_dhat(), params.r);
    double a4 =
        PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, double>::mollifier(
            x, dist * dist);

    return a1 * a2 * a3 * a4;
}

template <int max_vert, typename PrimitiveA, typename PrimitiveB>
Vector<
    double,
    -1,
    SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>
SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::gradient(
    const Vector<
        double,
        -1,
        SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>&
        positions,
    const ParameterType& params) const
{
    Eigen::VectorXi core_indices(n_core_dofs);
    for (int i = 0; i < n_core_dofs_A; i++)
        core_indices(i) = i;
    for (int i = 0; i < n_core_dofs_B; i++)
        core_indices(i + n_core_dofs_A) = i + pA->n_dofs();

    Vector<double, n_core_dofs> x;
    x = positions(core_indices);

    const auto dtype = PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_distance_type(x);
    
    Vector<double, dim> closest_direction;
    Eigen::Matrix<double, dim, n_core_dofs> closest_direction_grad;
    std::tie(closest_direction, closest_direction_grad) = PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_closest_direction_gradient(x, dtype);
    
    const double dist = closest_direction.norm();

    // these two use autodiff with different variable count
    auto gA_reduced = pA->grad(closest_direction, positions.head(pA->n_dofs()));
    auto gB_reduced = pB->grad(-closest_direction, positions.tail(pB->n_dofs()));

    // gradient of barrier potential
    double barrier = 0;
    Vector<double, n_core_dofs> gBarrier = Vector<double, n_core_dofs>::Zero();
    {
        barrier = Math<double>::inv_barrier(dist / Super::get_dhat(), params.r);

        const Vector<double, dim> closest_direction_normalized =
            closest_direction / dist;
        const double barrier_1st_deriv =
            Math<double>::inv_barrier_grad(dist / Super::get_dhat(), params.r)
            / Super::get_dhat();
        const Vector<double, dim> gBarrier_wrt_d =
            barrier_1st_deriv * closest_direction_normalized;
        gBarrier = closest_direction_grad.transpose() * gBarrier_wrt_d;
    }

    // gradient of mollifier
    {
        double mollifier = 0;
        Vector<double, n_core_dofs> gMollifier = Vector<double, n_core_dofs>::Zero();
#ifdef DERIVATIVES_WITH_AUTODIFF
        DiffScalarBase::setVariableCount(n_core_dofs);
        using T = ADGrad<n_core_dofs>;
        Vector<T, n_core_dofs> xAD = slice_positions<T, n_core_dofs, 1>(x);
        Vector<T, dim> closest_direction_autodiff = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::compute_closest_direction(xAD, dtype);
        const auto dist_sqr_AD = closest_direction_autodiff.squaredNorm();
        auto mollifier_autodiff =
            PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::mollifier(
                xAD, dist_sqr_AD);
        mollifier = mollifier_autodiff.getValue();
        gMollifier = mollifier_autodiff.getGradient();
#else
        Vector<double, n_core_dofs + 1> mollifier_grad;
        std::tie(mollifier, mollifier_grad) =
            PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_mollifier_gradient(x, dist * dist);

        const Vector<double, n_core_dofs> dist_sqr_grad = 2 * closest_direction_grad.transpose() * closest_direction;
        mollifier_grad.head(n_core_dofs) += mollifier_grad(n_core_dofs) * dist_sqr_grad;
        gMollifier = mollifier_grad.head(n_core_dofs);
#endif
        // merge mollifier into barrier
        gBarrier = gBarrier * mollifier + gMollifier * barrier;
        barrier *= mollifier;
    }

    // grad of tangent/normal terms
    double orient = 0;
    Vector<double, -1, max_size> gOrient;
    {
        Vector<double, -1, max_size> gA = Vector<double, -1, max_size>::Zero(ndofs()), gB = Vector<double, -1, max_size>::Zero(ndofs());
        {
            gA(core_indices) = closest_direction_grad.transpose() * gA_reduced.head(dim);
            gA.head(pA->n_dofs()) += gA_reduced.tail(pA->n_dofs());

            gB(core_indices) = closest_direction_grad.transpose() * -gB_reduced.head(dim);
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

template <int max_vert, typename PrimitiveA, typename PrimitiveB>
MatrixMax<
    double,
    SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size,
    SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>
SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::hessian(
    const Vector<
        double,
        -1,
        SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>&
        positions,
    const ParameterType& params) const
{
    Eigen::VectorXi core_indices(n_core_dofs);
    for (int i = 0; i < n_core_dofs_A; i++)
        core_indices(i) = i;
    for (int i = 0; i < n_core_dofs_B; i++)
        core_indices(i + n_core_dofs_A) = i + pA->n_dofs();

    Vector<double, n_core_dofs> x;
    x = positions(core_indices);

    const auto dtype = PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_distance_type(x);
    
    Vector<double, dim> closest_direction;
    Eigen::Matrix<double, dim, n_core_dofs> closest_direction_grad;
    std::array<Eigen::Matrix<double, n_core_dofs, n_core_dofs>, dim> closest_direction_hess;
    std::tie(closest_direction, closest_direction_grad, closest_direction_hess) = PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_closest_direction_hessian(x, dtype);
    
    const double dist = closest_direction.norm();

    // these two use autodiff with different variable count
    auto gA_reduced = pA->grad(closest_direction, positions.head(pA->n_dofs()));
    auto hA_reduced = pA->hessian(closest_direction, positions.head(pA->n_dofs()));
    auto gB_reduced = pB->grad(-closest_direction, positions.tail(pB->n_dofs()));
    auto hB_reduced = pB->hessian(-closest_direction, positions.tail(pB->n_dofs()));

    // hessian of barrier potential
    double barrier = 0;
    Vector<double, n_core_dofs> gBarrier = Vector<double, n_core_dofs>::Zero();
    Eigen::Matrix<double, n_core_dofs, n_core_dofs> hBarrier = Eigen::Matrix<double, n_core_dofs, n_core_dofs>::Zero();
    {
        barrier = Math<double>::inv_barrier(dist / Super::get_dhat(), params.r);

        const Vector<double, dim> closest_direction_normalized =
            closest_direction / dist;
        const double barrier_1st_deriv =
            Math<double>::inv_barrier_grad(dist / Super::get_dhat(), params.r)
            / Super::get_dhat();
        const Vector<double, dim> gBarrier_wrt_d =
            barrier_1st_deriv * closest_direction_normalized;
        gBarrier = closest_direction_grad.transpose() * gBarrier_wrt_d;

        const double barrier_2nd_deriv =
            Math<double>::inv_barrier_hess(dist / Super::get_dhat(), params.r)
            / Super::get_dhat() / Super::get_dhat();
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
        Vector<double, n_core_dofs> gMollifier = Vector<double, n_core_dofs>::Zero();
        Eigen::Matrix<double, n_core_dofs, n_core_dofs> hMollifier = Eigen::Matrix<double, n_core_dofs, n_core_dofs>::Zero();
#ifdef DERIVATIVES_WITH_AUTODIFF
        DiffScalarBase::setVariableCount(n_core_dofs);
        using T = ADHessian<n_core_dofs>;
        Vector<T, n_core_dofs> xAD = slice_positions<T, n_core_dofs, 1>(x);
        Vector<T, dim> closest_direction_autodiff = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::compute_closest_direction(xAD, dtype);
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
        std::tie(mollifier, mollifier_grad, mollifier_hess) =
            PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_mollifier_hessian(x, dist * dist);

        const Vector<double, n_core_dofs> dist_sqr_grad = 2 * closest_direction_grad.transpose() * closest_direction;
        mollifier_grad.head(n_core_dofs) += mollifier_grad(n_core_dofs) * dist_sqr_grad;
        Eigen::Matrix<double, n_core_dofs, n_core_dofs> dist_sqr_hess = 2 * closest_direction_grad.transpose() * closest_direction_grad;
        for (int d = 0; d < dim; d++)
            dist_sqr_hess += 2 * closest_direction(d) * closest_direction_hess[d];
        mollifier_hess.topLeftCorner(n_core_dofs, n_core_dofs) +=
            dist_sqr_hess * mollifier_grad(n_core_dofs)
            + dist_sqr_grad * mollifier_hess(n_core_dofs, n_core_dofs) * dist_sqr_grad.transpose()
            + dist_sqr_grad * mollifier_hess.block(n_core_dofs, 0, 1, n_core_dofs)
            + mollifier_hess.block(0, n_core_dofs, n_core_dofs, 1) * dist_sqr_grad.transpose();

        gMollifier = mollifier_grad.head(core_indices.size());
        hMollifier = mollifier_hess.topLeftCorner(core_indices.size(), core_indices.size());
#endif
        // merge mollifier into barrier
        hBarrier = hBarrier * mollifier + gBarrier * gMollifier.transpose() + gMollifier * gBarrier.transpose() + hMollifier * barrier;
        gBarrier = gBarrier * mollifier + gMollifier * barrier;
        barrier *= mollifier;
    }

    // grad of tangent/normal terms
    double orient = 0;
    Vector<double, -1, max_size> gOrient;
    MatrixMax<double, max_size, max_size> hOrient;
    {
        Vector<double, -1, max_size> gA = Vector<double, -1, max_size>::Zero(ndofs()), gB = Vector<double, -1, max_size>::Zero(ndofs());
        MatrixMax<double, max_size, max_size> hA = MatrixMax<double, max_size, max_size>::Zero(ndofs(), ndofs()), hB = MatrixMax<double, max_size, max_size>::Zero(ndofs(), ndofs());
        {
            gA(core_indices) = closest_direction_grad.transpose() * gA_reduced.head(dim);
            gA.head(pA->n_dofs()) += gA_reduced.tail(pA->n_dofs());

            hA(core_indices, core_indices) = closest_direction_grad.transpose() * hA_reduced.topLeftCorner(dim, dim) * closest_direction_grad;
            for (int d = 0; d < dim; d++)
                hA(core_indices, core_indices) += gA_reduced(d) * closest_direction_hess[d];

            hA.topLeftCorner(pA->n_dofs(), pA->n_dofs()) +=
                hA_reduced.bottomRightCorner(pA->n_dofs(), pA->n_dofs());
            hA(core_indices, Eigen::seqN(0, pA->n_dofs())) += closest_direction_grad.transpose()
                * hA_reduced.topRightCorner(dim, pA->n_dofs());
            hA(Eigen::seqN(0, pA->n_dofs()), core_indices) +=
                hA_reduced.bottomLeftCorner(pA->n_dofs(), dim)
                * closest_direction_grad;

            gB(core_indices) = closest_direction_grad.transpose() * -gB_reduced.head(dim);
            gB.tail(pB->n_dofs()) += gB_reduced.tail(pB->n_dofs());

            hB(core_indices, core_indices) = closest_direction_grad.transpose() * hB_reduced.topLeftCorner(dim, dim)
                * closest_direction_grad;
            for (int d = 0; d < dim; d++)
                hB(core_indices, core_indices) -= gB_reduced(d) * closest_direction_hess[d];

            hB.bottomRightCorner(pB->n_dofs(), pB->n_dofs()) +=
                hB_reduced.bottomRightCorner(pB->n_dofs(), pB->n_dofs());
            hB(core_indices, Eigen::seqN(pA->n_dofs(), pB->n_dofs())) -= closest_direction_grad.transpose()
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

    return hOrient;
}

// ---- distance ----

template <int max_vert, typename PrimitiveA, typename PrimitiveB>
double
SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::compute_distance(
    const Vector<
        double,
        -1,
        SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>&
        positions) const
{
    Vector<double, n_core_points * dim> x;
    x << positions.head(PrimitiveA::n_core_points * dim),
        positions.segment(pA->n_dofs(), PrimitiveB::n_core_points * dim);

    // grad of "d" wrt. points
    Vector<double, dim> closest_direction =
        PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, double>::
            compute_closest_direction(x, DTYPE::AUTO);

    return closest_direction.squaredNorm();
}

template <int max_vert, typename PrimitiveA, typename PrimitiveB>
Vector<
    double,
    -1,
    SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>
SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::
    compute_distance_gradient(
        const Vector<
            double,
            -1,
            SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::
                max_size>& positions) const
{
    return Vector<
        double, -1,
        SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>::
        Zero(ndofs());
}

template <int max_vert, typename PrimitiveA, typename PrimitiveB>
MatrixMax<
    double,
    SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size,
    SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>
SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::
    compute_distance_hessian(
        const Vector<
            double,
            -1,
            SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::
                max_size>& positions) const
{
    return MatrixMax<
        double,
        SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size,
        SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>::
        Zero(ndofs(), ndofs());
}

// Note: Primitive pair order cannot change
template class SmoothCollisionTemplate<max_vert_3d, Edge3, Point3>;
template class SmoothCollisionTemplate<max_vert_3d, Edge3, Edge3>;
template class SmoothCollisionTemplate<max_vert_3d, Point3, Point3>;
template class SmoothCollisionTemplate<max_vert_3d, Face, Point3>;
} // namespace ipc