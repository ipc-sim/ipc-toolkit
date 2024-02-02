#include "smooth_collision.hpp"
#include <ipc/smooth_contact/primitives/point.hpp>
#include <ipc/smooth_contact/primitives/edge.hpp>
#include <ipc/smooth_contact/primitives/face.hpp>

namespace ipc {
    namespace {
        template <int max_vert, typename T, int dim, int max_dim = dim>
        Vector<T, dim, max_dim> reorder_vector(const SmoothCollision<max_vert> &a, const SmoothCollision<max_vert> &b, const Vector<T, dim, max_dim> &vec_a)
        {
            auto ids_a = a.vertices;
            auto ids_b = b.vertices;

            auto vec_b = vec_a;
            for (int i = 0; i < ids_a.size(); i++)
            {
                if (ids_a[i] < 0)
                    break;
                for (int j = 0; j < ids_b.size(); j++)
                {
                    if (ids_b[j] == ids_a[i])
                    {
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
        MatrixMax<T, dim, dim> reorder_matrix(const SmoothCollision<max_vert> &a, const SmoothCollision<max_vert> &b, const MatrixMax<T, dim, dim>&mat_a)
        {
            auto ids_a = a.vertices;
            auto ids_b = b.vertices;

            std::map<int, int> id_map;
            for (int i = 0; i < ids_a.size(); i++)
            {
                if (ids_a[i] < 0)
                    break;
                for (int j = 0; j < ids_b.size(); j++)
                {
                    if (ids_b[j] == ids_a[i])
                    {
                        id_map[i] = j;
                        break;
                    }
                    if (ids_b[j] < 0)
                        break;
                }
            }

            auto mat_b = mat_a;
            for (int i = 0; i < ids_a.size(); i++)
                for (int j = 0; j < ids_a.size(); j++)
                {
                    if (ids_a[i] < 0 || ids_a[j] < 0)
                        break;
                    mat_b.block(id_map[i] * 3, id_map[j] * 3, 3, 3) = mat_a.block(i * 3, j * 3, 3, 3);
                }

            return mat_b;
        }
    }

    template <int max_vert, typename PrimitiveA, typename PrimitiveB>
    SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::SmoothCollisionTemplate(
        long primitive0_,
        long primitive1_,
        SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::DTYPE dtype,
        const CollisionMesh &mesh,
        const ParameterType &param,
        const double &dhat,
        const Eigen::MatrixXd &V)
    : SmoothCollision<max_vert>(primitive0_, primitive1_, dhat, mesh)
    {
        VectorMax3d d = PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_closest_direction(mesh, V, primitive0_, primitive1_, dtype);
        pA = std::make_shared<PrimitiveA>(primitive0_, mesh, V,  d, param.alpha, param.beta);
        pB = std::make_shared<PrimitiveB>(primitive1_, mesh, V, -d, param.alpha, param.beta);
        int i = 0;
        for (auto& v : pA->vertex_ids())
            Super::vertices[i++] = v;
        for (auto& v : pB->vertex_ids())
            Super::vertices[i++] = v;
        assert(i == pA->n_vertices() + pB->n_vertices());
        Super::is_active_ = (d.norm() < Super::get_dhat()) && pA->is_active() && pB->is_active();

        if (!Super::is_active())
            return;

        // test derivatives
        // std::shared_ptr<SmoothCollision<max_vert>> debug_ptr;
        // if constexpr (std::is_same<PrimitiveA, Edge3>::value && std::is_same<PrimitiveB, Edge3>::value)
        //     debug_ptr = std::make_shared<SmoothEdgeEdge3Collision>(Super::primitive0, Super::primitive1, mesh, param, dhat, V);
        // else if constexpr (std::is_same<PrimitiveA, Face>::value && std::is_same<PrimitiveB, Point3>::value)
        //     debug_ptr = std::make_shared<SmoothFaceVertexCollision>(Super::primitive0, Super::primitive1, mesh, param, dhat, V);
        // else if constexpr (std::is_same<PrimitiveA, Edge3>::value && std::is_same<PrimitiveB, Point3>::value)
        //     debug_ptr = std::make_shared<SmoothEdgeVertex3Collision>(Super::primitive0, Super::primitive1, mesh, param, dhat, V);
        // else if constexpr (std::is_same<PrimitiveA, Point3>::value && std::is_same<PrimitiveB, Point3>::value)
        //     debug_ptr = std::make_shared<SmoothVertexVertex3Collision>(Super::primitive0, Super::primitive1, mesh, param, dhat, V);
        
        // if (!debug_ptr)
        //     return;

        // double val = (*this)(this->dof(V, mesh.edges(), mesh.faces()), param);
        // double cc_val = (*debug_ptr)(debug_ptr->dof(V, mesh.edges(), mesh.faces()), param);
        
        // if (std::abs(cc_val - val) > 1e-6 * std::min(std::abs(cc_val), std::abs(val)))
        //     logger().error("Inconsistent potential: {} vs {}", cc_val, val);
    }

    template <int max_vert, typename PrimitiveA, typename PrimitiveB>
    SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::~SmoothCollisionTemplate()
    {
    }

    template <int max_vert, typename PrimitiveA, typename PrimitiveB>
    double SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::operator()(
        const Vector<double, -1, SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>& positions, 
        const ParameterType &params) const
    {
        Vector<double, n_core_points*dim> x;
        x << positions.head(PrimitiveA::n_core_points*dim), positions.segment(pA->n_dofs(), PrimitiveB::n_core_points*dim);
        
        // grad of "d" wrt. points
        const Vector<double, dim> closest_direction = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, double>::compute_closest_direction(x, DTYPE::AUTO);
        const double dist = closest_direction.norm();

        assert(positions.size() == pA->n_dofs() + pB->n_dofs());
        double a1 = pA->potential(closest_direction, positions.head(pA->n_dofs()));
        double a2 = pB->potential(-closest_direction, positions.tail(pB->n_dofs()));
        double a3 = Math<double>::inv_barrier(dist / Super::get_dhat(), params.r);
        double a4 = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, double>::mollifier(x, dist * dist);

        return a1 * a2 * a3 * a4;
    }

    template <int max_vert, typename PrimitiveA, typename PrimitiveB>
    Vector<double, -1, SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size> SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::gradient(
        const Vector<double, -1, SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>& positions, 
        const ParameterType &params) const
    {
        Vector<double, n_core_dofs> x_double;
        x_double << positions.head(n_core_dofs_A), positions.segment(pA->n_dofs(), n_core_dofs_B);

        const auto dtype = PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_distance_type(x_double);
        const Vector<double, dim> closest_direction = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, double>::compute_closest_direction(x_double, dtype);
        const double dist = closest_direction.norm();

        // these two use autodiff with different variable count
        auto gA_reduced = pA->grad(closest_direction, positions.head(pA->n_dofs()));
        auto gB_reduced = pB->grad(-closest_direction, positions.tail(pB->n_dofs()));

        DiffScalarBase::setVariableCount(n_core_dofs);
        using T = ADGrad<n_core_dofs>;

        Vector<T, n_core_dofs> x = slice_positions<T, n_core_dofs, 1>(x_double);

        MatrixMax<double, max_size, dim> closest_direction_grad;
        closest_direction_grad.setZero(ndofs(), dim);
        auto closest_direction_autodiff = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::compute_closest_direction(x, dtype);
        for (int k = 0; k < dim; k++)
        {
            closest_direction_grad.block(0, k, n_core_dofs_A, 1) = closest_direction_autodiff(k).getGradient().head(n_core_dofs_A);
            closest_direction_grad.block(pA->n_dofs(), k, n_core_dofs_B, 1) = closest_direction_autodiff(k).getGradient().tail(n_core_dofs_B);
        }

        double out = 0.;
        Vector<double, -1, max_size> gOut(ndofs());
        gOut.setZero();

        // grad of tangent/normal terms
        {
            const double potential_a = pA->potential(closest_direction, positions.head(pA->n_dofs()));
            const double potential_b = pB->potential(-closest_direction, positions.tail(pB->n_dofs()));

            gOut += closest_direction_grad * (gA_reduced.head(dim) * potential_b);
            gOut.head(pA->n_dofs()) += gA_reduced.tail(pA->n_dofs()) * potential_b;
            
            gOut += closest_direction_grad * (-gB_reduced.head(dim) * potential_a);
            gOut.tail(pB->n_dofs()) += gB_reduced.tail(pB->n_dofs()) * potential_a;
            
            out = potential_a * potential_b;
        }

        // grad of barrier potential
        {
            double barrier_deriv = Math<double>::inv_barrier_grad(dist / Super::get_dhat(), params.r);
            Vector<double, -1, max_size> gBarrier = (barrier_deriv / Super::get_dhat() / dist) * closest_direction;
            double barrier = Math<double>::inv_barrier(dist / Super::get_dhat(), params.r);

            gOut = gOut * barrier + closest_direction_grad * gBarrier * out;
            out *= barrier;
        }

        // grad of mollifier
        {
            auto mollifier_autodiff = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::mollifier(x, closest_direction_autodiff.squaredNorm());
            
            gOut *= mollifier_autodiff.getValue();
            gOut.head(n_core_dofs_A) += mollifier_autodiff.getGradient().head(n_core_dofs_A) * out;
            gOut.segment(pA->n_dofs(), n_core_dofs_B) += mollifier_autodiff.getGradient().tail(n_core_dofs_B) * out;
            out *= mollifier_autodiff.getValue();
        }

        return gOut;
    }

    template <int max_vert, typename PrimitiveA, typename PrimitiveB>
    MatrixMax<double, SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size, SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size> SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::hessian(
        const Vector<double, -1, SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>& positions, 
        const ParameterType &params) const
    {
        Vector<double, n_core_dofs> x_double;
        x_double << positions.head(n_core_dofs_A), positions.segment(pA->n_dofs(), n_core_dofs_B);

        const auto dtype = PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_distance_type(x_double);
        const Vector<double, dim> closest_direction = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, double>::compute_closest_direction(x_double, dtype);
        const double dist = closest_direction.norm();

        // these two use autodiff with different variable count
        auto gA_reduced = pA->grad(closest_direction, positions.head(pA->n_dofs()));
        auto hA_reduced = pA->hessian(closest_direction, positions.head(pA->n_dofs()));
        auto gB_reduced = pB->grad(-closest_direction, positions.tail(pB->n_dofs()));
        auto hB_reduced = pB->hessian(-closest_direction, positions.tail(pB->n_dofs()));

        DiffScalarBase::setVariableCount(n_core_dofs);
        using T = ADHessian<n_core_dofs>;

        Vector<T, n_core_dofs> x = slice_positions<T, n_core_dofs, 1>(x_double);

        auto closest_direction_autodiff = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::compute_closest_direction(x, dtype);
        Eigen::Matrix<double, -1, dim, Eigen::ColMajor, max_size, dim> closest_direction_grad;
        closest_direction_grad.setZero(ndofs(), dim);
        std::array<MatrixMax<double, max_size, max_size>, dim> closest_direction_hess;
        for (int k = 0; k < dim; k++)
        {
            closest_direction_grad.block(0, k, n_core_dofs_A, 1) = closest_direction_autodiff(k).getGradient().head(n_core_dofs_A);
            closest_direction_grad.block(pA->n_dofs(), k, n_core_dofs_B, 1) = closest_direction_autodiff(k).getGradient().tail(n_core_dofs_B);

            closest_direction_hess[k].setZero(ndofs(), ndofs());
            const auto& hess = closest_direction_autodiff(k).getHessian();
            closest_direction_hess[k].block(0, 0, n_core_dofs_A, n_core_dofs_A) = hess.topLeftCorner(n_core_dofs_A, n_core_dofs_A);
            closest_direction_hess[k].block(pA->n_dofs(), pA->n_dofs(), n_core_dofs_B, n_core_dofs_B) = hess.bottomRightCorner(n_core_dofs_B, n_core_dofs_B);
            closest_direction_hess[k].block(0, pA->n_dofs(), n_core_dofs_A, n_core_dofs_B) = hess.topRightCorner(n_core_dofs_A, n_core_dofs_B);
            closest_direction_hess[k].block(pA->n_dofs(), 0, n_core_dofs_B, n_core_dofs_A) = hess.bottomLeftCorner(n_core_dofs_B, n_core_dofs_A);
        }

        // grad of tangent/normal terms
        double orient = 0;
        Vector<double, -1, max_size> gOrient;
        MatrixMax<double, max_size, max_size> hOrient;
        {
            Vector<double, -1, max_size> gA, gB;
            MatrixMax<double, max_size, max_size> hA, hB;
            {
                gA = closest_direction_grad * gA_reduced.head(dim);
                gA.head(pA->n_dofs()) += gA_reduced.tail(pA->n_dofs());

                hA = closest_direction_grad * hA_reduced.topLeftCorner(dim, dim) * closest_direction_grad.transpose();
                for (int d = 0; d < dim; d++)
                    hA += gA_reduced(d) * closest_direction_hess[d];

                hA.topLeftCorner(pA->n_dofs(), pA->n_dofs()) += hA_reduced.bottomRightCorner(pA->n_dofs(), pA->n_dofs());
                hA.leftCols(pA->n_dofs()) += closest_direction_grad * hA_reduced.topRightCorner(dim, pA->n_dofs());
                hA.topRows(pA->n_dofs()) += hA_reduced.bottomLeftCorner(pA->n_dofs(), dim) * closest_direction_grad.transpose();

                gB = closest_direction_grad * -gB_reduced.head(dim);
                gB.tail(pB->n_dofs()) += gB_reduced.tail(pB->n_dofs());

                hB = closest_direction_grad * hB_reduced.topLeftCorner(dim, dim) * closest_direction_grad.transpose();
                for (int d = 0; d < dim; d++)
                    hB -= gB_reduced(d) * closest_direction_hess[d];

                hB.bottomRightCorner(pB->n_dofs(), pB->n_dofs()) += hB_reduced.bottomRightCorner(pB->n_dofs(), pB->n_dofs());
                hB.rightCols(pB->n_dofs()) -= closest_direction_grad * hB_reduced.topRightCorner(dim, pB->n_dofs());
                hB.bottomRows(pB->n_dofs()) -= hB_reduced.bottomLeftCorner(pB->n_dofs(), dim) * closest_direction_grad.transpose();
            }
            const double potential_a = pA->potential(closest_direction, positions.head(pA->n_dofs()));
            const double potential_b = pB->potential(-closest_direction, positions.tail(pB->n_dofs()));
            
            orient = potential_a * potential_b;
            gOrient = gA * potential_b + gB * potential_a;
            hOrient = (potential_a * hB + potential_b * hA) + (gA * gB.transpose() + gB * gA.transpose());
        }

        // grad of barrier potential
        double barrier = 0;
        Vector<double, -1, max_size> gBarrier;
        MatrixMax<double, max_size, max_size> hBarrier;
        {
            barrier = Math<double>::inv_barrier(dist / Super::get_dhat(), params.r);

            const Vector<double, dim> closest_direction_normalized = closest_direction / dist;
            const double barrier_1st_deriv = Math<double>::inv_barrier_grad(dist / Super::get_dhat(), params.r) / Super::get_dhat();
            const Vector<double, dim> gBarrier_wrt_d = barrier_1st_deriv * closest_direction_normalized;
            gBarrier = closest_direction_grad * gBarrier_wrt_d;

            const double barrier_2nd_deriv = Math<double>::inv_barrier_hess(dist / Super::get_dhat(), params.r) / Super::get_dhat() / Super::get_dhat();
            const Eigen::Matrix<double, dim, dim> hBarrier_wrt_d = (barrier_1st_deriv / dist) * Eigen::Matrix<double, dim, dim>::Identity() + 
                    (barrier_2nd_deriv - barrier_1st_deriv / dist) * closest_direction_normalized * closest_direction_normalized.transpose();
            hBarrier = closest_direction_grad * hBarrier_wrt_d * closest_direction_grad.transpose();
            for (int d = 0; d < dim; d++)
                hBarrier += closest_direction_hess[d] * gBarrier_wrt_d(d);
        }

        // grad of mollifier
        double mollifier = 0;
        Vector<double, -1, max_size> gMollifier;
        gMollifier.setZero(ndofs());
        MatrixMax<double, max_size, max_size> hMollifier;
        hMollifier.setZero(ndofs(), ndofs());
        {
            auto mollifier_autodiff = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::mollifier(x, closest_direction_autodiff.squaredNorm());
            mollifier = mollifier_autodiff.getValue();
            gMollifier.segment(0, n_core_dofs_A) = mollifier_autodiff.getGradient().head(n_core_dofs_A);
            gMollifier.segment(pA->n_dofs(), n_core_dofs_B) = mollifier_autodiff.getGradient().tail(n_core_dofs_B);

            hMollifier.block(0, 0, n_core_dofs_A, n_core_dofs_A) = mollifier_autodiff.getHessian().topLeftCorner(n_core_dofs_A, n_core_dofs_A);
            hMollifier.block(pA->n_dofs(), pA->n_dofs(), n_core_dofs_B, n_core_dofs_B) = mollifier_autodiff.getHessian().bottomRightCorner(n_core_dofs_B, n_core_dofs_B);
            hMollifier.block(0, pA->n_dofs(), n_core_dofs_A, n_core_dofs_B) = mollifier_autodiff.getHessian().topRightCorner(n_core_dofs_A, n_core_dofs_B);
            hMollifier.block(pA->n_dofs(), 0, n_core_dofs_B, n_core_dofs_A) = mollifier_autodiff.getHessian().bottomLeftCorner(n_core_dofs_B, n_core_dofs_A);
        }

        DiffScalarBase::setVariableCount(ndofs());
        ADHessian<-1, max_size> orient_term(orient, gOrient, hOrient);
        ADHessian<-1, max_size> mollifier_term(mollifier, gMollifier, hMollifier);
        ADHessian<-1, max_size> barrier_term(barrier, gBarrier, hBarrier);
        MatrixMax<double, max_size, max_size> out = (orient_term * mollifier_term * barrier_term).getHessian();

        return out;
    }

    // ---- distance ----

    template <int max_vert, typename PrimitiveA, typename PrimitiveB> double 
    SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::compute_distance(const Vector<double, -1, SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>& positions) const
    {
        Vector<double, n_core_points*dim> x;
        x << positions.head(PrimitiveA::n_core_points*dim), positions.segment(pA->n_dofs(), PrimitiveB::n_core_points*dim);
        
        // grad of "d" wrt. points
        Vector<double, dim> closest_direction = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, double>::compute_closest_direction(x, DTYPE::AUTO);

        return closest_direction.squaredNorm();
    }

    template <int max_vert, typename PrimitiveA, typename PrimitiveB> Vector<double, -1, SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>
    SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::compute_distance_gradient(const Vector<double, -1, SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>& positions) const
    {
        return Vector<double, -1, SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>::Zero(ndofs());
    }

    template <int max_vert, typename PrimitiveA, typename PrimitiveB> MatrixMax<double, SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size, SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>
    SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::compute_distance_hessian(const Vector<double, -1, SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>& positions) const
    {
        return MatrixMax<double, SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size, SmoothCollisionTemplate<max_vert, PrimitiveA, PrimitiveB>::max_size>::Zero(ndofs(), ndofs());
    }

    // Note: Primitive pair order cannot change
    template class SmoothCollisionTemplate<max_vert_3d, Edge3 , Point3>;
    template class SmoothCollisionTemplate<max_vert_3d, Edge3 , Edge3 >;
    template class SmoothCollisionTemplate<max_vert_3d, Point3, Point3>;
    template class SmoothCollisionTemplate<max_vert_3d, Face  , Point3>;
}