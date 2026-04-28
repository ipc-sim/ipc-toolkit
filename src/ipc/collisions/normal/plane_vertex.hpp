#pragma once

#include <ipc/candidates/plane_vertex.hpp>
#include <ipc/collisions/normal/normal_collision.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class PlaneVertexNormalCollision : public PlaneVertexCandidate,
                                   public NormalCollision {
public:
    using PlaneVertexCandidate::PlaneVertexCandidate;

    PlaneVertexNormalCollision(const PlaneVertexCandidate& candidate)
        : PlaneVertexCandidate(candidate)
    {
    }

    PlaneVertexNormalCollision(
        const Eigen::Hyperplane<double, 3>& _plane,
        const index_t _vertex_id,
        const double _weight,
        const Eigen::SparseVector<double>& _weight_gradient)
        : PlaneVertexCandidate(_plane, _vertex_id)
        , NormalCollision(_weight, _weight_gradient)
    {
    }
};

} // namespace ipc
