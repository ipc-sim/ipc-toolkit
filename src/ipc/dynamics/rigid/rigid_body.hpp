#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc::rigid {

struct RigidBody {
    VectorMax3d rotation_vector; // Affine matrix
    VectorMax3d position;        // Translation vector

    RigidBody() = default;

    RigidBody(const VectorMax3d& r, const VectorMax3d& p)
        : rotation_vector(r)
        , position(p)
    {
    }

    Eigen::MatrixXd
    transform_vertices(const Eigen::MatrixXd& rest_positions) const;
};

} // namespace ipc::rigid