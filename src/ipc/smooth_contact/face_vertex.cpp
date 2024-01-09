#include "face_vertex.hpp"
#include "smooth_point_face.hpp"
#include <ipc/utils/AutodiffTypes.hpp>
#include <ipc/utils/logger.hpp>

namespace ipc {
    namespace {
        template <class T, int dim = 3, int n_pts = 4>
        std::array<Vector<T, dim>, n_pts> slice_positions(const VectorMax12d &positions)
        {
            std::array<Vector<T, dim>, n_pts> out;
            for (int d = 0; d < dim; d++)
                for (int i = 0; i < n_pts; i++)
                    out[i](d) = T(i * dim + d, i * dim + positions(d));

            return out;
        }
    }

    double SmoothFaceVertexCollision::operator()(const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        return smooth_point_face_potential_single_point<double>(positions.head<3>(), positions.segment<3>(3), positions.segment<3>(6), positions.tail<3>(), params, known_dtype());
    }

    VectorMax12d SmoothFaceVertexCollision::gradient(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarGrad<12>;
        auto [p, v0, v1, v2] = slice_positions<Diff>(positions);

        const auto val = smooth_point_face_potential_single_point<Diff>(p, v0, v1, v2, params, known_dtype());

        return val.getGradient();
    }

    MatrixMax12d SmoothFaceVertexCollision::hessian(
        const VectorMax12d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarHessian<12>;
        auto [p, v0, v1, v2] = slice_positions<Diff>(positions);

        const auto val = smooth_point_face_potential_single_point<Diff>(p, v0, v1, v2, params, known_dtype());

        return val.getHessian();
    }
}