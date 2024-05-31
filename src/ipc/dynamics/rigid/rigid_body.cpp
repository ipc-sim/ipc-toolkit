#include "rigid_body.hpp"

#include <ipc/utils/sinc.hpp>

namespace ipc::rigid {

namespace {
    inline Eigen::Matrix3d cross_product_matrix(const Eigen::Vector3d& x)
    {
        Eigen::Matrix3d X;
        X << 0, -x.z(), x.y(), //
            x.z(), 0, -x.x(),  //
            -x.y(), x.x(), 0;
        return X;
    }
} // namespace

Eigen::MatrixXd
RigidBody::transform_vertices(const Eigen::MatrixXd& rest_positions) const
{
    // Compute: R(θ) x̄ + p
    // transpose because x is row-ordered
    const int dim = rest_positions.cols();

    assert(
        (dim == 2 && rotation_vector.size() == 1)
        || (dim == 3 && rotation_vector.size() == 3));

    // Convert the rotation vector to a rotation matrix
    MatrixMax3d R(dim, dim);
    if (dim == 2) {
        const double theta = rotation_vector(0);
        R << sin(theta), -cos(theta), cos(theta), sin(theta);
    } else {
        assert(dim == 3);
        const double sinc_angle = sinc_normx(rotation_vector);
        const double sinc_half_angle = sinc_normx((rotation_vector / 2).eval());
        const Eigen::Matrix3d K = cross_product_matrix(rotation_vector);
        const Eigen::Matrix3d K2 = K * K;
        R = sinc_angle * K + 0.5 * sinc_half_angle * sinc_half_angle * K2;
        R.diagonal().array() += 1.0;
    }

    return (rest_positions * R.transpose()).rowwise() + position.transpose();
}

} // namespace ipc::rigid