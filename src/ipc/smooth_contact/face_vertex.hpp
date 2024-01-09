#pragma once

#include <ipc/collisions/face_vertex.hpp>

namespace ipc {

class SmoothFaceVertexCollision : public FaceVertexCollision {
public:
    using FaceVertexCollision::FaceVertexCollision;

    double operator()(const VectorMax12d& positions, 
        const ParameterType &params) const override;

    VectorMax12d gradient(
        const VectorMax12d& positions, 
        const ParameterType &params) const override;

    MatrixMax12d hessian(
        const VectorMax12d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd = false) const override;
};

} // namespace ipc
