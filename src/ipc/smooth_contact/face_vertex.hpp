#pragma once

#include <ipc/collisions/face_vertex.hpp>

namespace ipc {

class SmoothFaceVertexCollision : public FaceVertexCollision {
public:
    using FaceVertexCollision::FaceVertexCollision;

    SmoothFaceVertexCollision(
        const long _face_id,
        const long _vertex_id,
        const double _weight,
        const Eigen::SparseVector<double>& _weight_gradient,
        const PointTriangleDistanceType &_dtype)
        : FaceVertexCollision(_face_id, _vertex_id, _weight, _weight_gradient)
        , dtype(_dtype)
    {
    }

    PointTriangleDistanceType known_dtype() const override
    {
        // The distance type is known because of Collisions::build()
        return dtype;
    }

    double operator()(const VectorMax12d& positions, 
        const ParameterType &params) const override;

    VectorMax12d gradient(
        const VectorMax12d& positions, 
        const ParameterType &params) const override;

    MatrixMax12d hessian(
        const VectorMax12d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd = false) const override;
    
private:
    mutable PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO;
};

} // namespace ipc
