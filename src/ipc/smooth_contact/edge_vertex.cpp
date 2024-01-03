#include "edge_vertex.hpp"
#include "smooth_point_edge.hpp"
#include <ipc/utils/AutodiffTypes.hpp>
DECLARE_DIFFSCALAR_BASE();

namespace ipc {
    namespace {
        template <class T>
        std::array<VectorMax3<T>, 3> slice_positions(const VectorMax12d &positions, const int &dim)
        {
            VectorMax3<T> p, e0, e1;
            p.setZero(dim);
            e0.setZero(dim);
            e1.setZero(dim);
            for (int d = 0; d < dim; d++)
            {
                p(d) = T(d, positions(d));
                e0(d) = T(dim + d, positions(dim + d));
                e1(d) = T(2*dim + d, positions(2*dim + d));
            }

            return {{p, e0, e1}};
        }
    }

    double SmoothEdgeVertexCollision::operator()(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        const int dim = positions.size() / num_vertices();
        assert(dim * num_vertices() == positions.size());
        return smooth_point_edge_potential_single_point<double>(positions.segment(0, dim), positions.segment(dim, dim), positions.segment(dim * 2, dim), params);
    }

    VectorMax12d SmoothEdgeVertexCollision::gradient(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        const int dim = positions.size() / num_vertices();
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarGrad<12>;
        auto [p, e0, e1] = slice_positions<Diff>(positions, dim);

        const auto val = smooth_point_edge_potential_single_point<Diff>(p, e0, e1, params);

        VectorMax12d grad;
        grad = val.getGradient().head(3*dim);

        return grad;
    }

    MatrixMax12d SmoothEdgeVertexCollision::hessian(
        const VectorMax12d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        const int dim = positions.size() / num_vertices();
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarHessian<12>;
        auto [p, e0, e1] = slice_positions<Diff>(positions, dim);

        const auto val = smooth_point_edge_potential_single_point<Diff>(p, e0, e1, params);

        MatrixMax12d hess;
        hess = val.getHessian().topLeftCorner(3*dim, 3*dim);

        return hess;
    }
}