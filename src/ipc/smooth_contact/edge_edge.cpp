#include "edge_edge.hpp"
#include "smooth_point_edge.hpp"
#include <ipc/utils/AutodiffTypes.hpp>

namespace ipc {
    double SmoothEdgeEdgeCollision::operator()(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        return 0.;
    }

    VectorMax12d SmoothEdgeEdgeCollision::gradient(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        VectorMax12d grad;
        grad.setZero(positions.size());
        return grad;
    }

    MatrixMax12d SmoothEdgeEdgeCollision::hessian(
        const VectorMax12d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {   
        MatrixMax12d hess;
        hess.setZero(positions.size(), positions.size());
        return hess;
    }
}