#include "edge_vertex.hpp"
#include <ipc/smooth_contact/smooth_point_edge.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
DECLARE_DIFFSCALAR_BASE();

namespace ipc {
    namespace {
        template <class T>
        std::tuple<VectorMax3<T>, VectorMax3<T>, VectorMax3<T>> slice_positions(const VectorMax12d &positions, const int &dim)
        {
            VectorMax3<T> p, e0, e1;
            p.setZero(dim);
            e0.setZero(dim);
            e1.setZero(dim);
            for (int d = 0; d < dim; d++)
            {
                p(d) = T(d, positions(d));
                e0(d) = T(3 + d, positions(dim + d));
                e1(d) = T(6 + d, positions(2*dim + d));
            }

            return std::make_tuple(p, e0, e1);
        }
    }

    double EdgeVertexCollision::operator()(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        const int dim = positions.size() / num_vertices();
        assert(dim * num_vertices() == positions.size());
        return smooth_point_edge_potential_single_point<double>(positions.segment(0, dim), positions.segment(dim, dim), positions.segment(dim * 2, dim), params.eps, params.alpha, params.r);
    }

    VectorMax12d EdgeVertexCollision::gradient(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        const int dim = positions.size() / num_vertices();
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarGrad<12>;
        auto [p, e0, e1] = slice_positions<Diff>(positions, dim);

        const auto val = smooth_point_edge_potential_single_point<Diff>(p, e0, e1, params.eps, params.alpha, params.r);

        VectorMax12d grad;
        grad.setZero(3 * dim);
        for (int d = 0; d < dim; d++)
            for (int i = 0; i < 3; i++)
                grad(d + dim * i) = val.getGradient()(d + 3 * i);
        
        return grad;
    }

    MatrixMax12d EdgeVertexCollision::hessian(
        const VectorMax12d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        const int dim = positions.size() / num_vertices();
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarHessian<12>;
        auto [p, e0, e1] = slice_positions<Diff>(positions, dim);

        const auto val = smooth_point_edge_potential_single_point<Diff>(p, e0, e1, params.eps, params.alpha, params.r);

        MatrixMax12d hess;
        hess.setZero(3 * dim, 3 * dim);
        for (int d1 = 0; d1 < dim; d1++)
            for (int i1 = 0; i1 < 3; i1++)
                for (int d2 = 0; d2 < dim; d2++)
                    for (int i2 = 0; i2 < 3; i2++)
                        hess(d1 + dim * i1, d2 + dim * i2) = val.getHessian()(d1 + 3 * i1, d2 + 3 * i2);
        
        return hess;
    }
}