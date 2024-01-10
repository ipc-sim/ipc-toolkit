#include "edge_vertex.hpp"
#include "smooth_point_edge.hpp"
#include <ipc/utils/AutodiffTypes.hpp>

DECLARE_DIFFSCALAR_BASE();

namespace ipc {
    namespace {
        template <class T>
        std::array<Vector2<T>, 3> slice_positions(const VectorMax12d &positions)
        {
            Vector2<T> p, e0, e1;
            for (int d = 0; d < 2; d++)
            {
                p(d) = T(d, positions(d));
                e0(d) = T(2 + d, positions(2 + d));
                e1(d) = T(4 + d, positions(4 + d));
            }

            return {{p, e0, e1}};
        }
    }

    double SmoothEdgeVertexCollision::operator()(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        const ParameterType local_params(local_eps, params.alpha, params.a, params.r, params.n_quadrature);

        assert(positions.size() == 6);
        return smooth_point_edge_potential_single_point<double>(positions.segment<2>(0), positions.segment<2>(2), positions.segment<2>(4), local_params);
   }

    VectorMax12d SmoothEdgeVertexCollision::gradient(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        const ParameterType local_params(local_eps, params.alpha, params.a, params.r, params.n_quadrature);

        assert(positions.size() == 6);
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarGrad<12>;
        auto [p, e0, e1] = slice_positions<Diff>(positions);

        const auto val = smooth_point_edge_potential_single_point<Diff>(p, e0, e1, local_params);

        VectorMax12d grad;
        grad = val.getGradient().head(6);

        return grad;
    }

    MatrixMax12d SmoothEdgeVertexCollision::hessian(
        const VectorMax12d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        const ParameterType local_params(local_eps, params.alpha, params.a, params.r, params.n_quadrature);

        assert(positions.size() == 6);
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarHessian<12>;
        auto [p, e0, e1] = slice_positions<Diff>(positions);

        const auto val = smooth_point_edge_potential_single_point<Diff>(p, e0, e1, local_params);

        MatrixMax12d hess;
        hess = val.getHessian().topLeftCorner(6, 6);

        return hess;
    }
}