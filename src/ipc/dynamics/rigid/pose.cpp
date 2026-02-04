#include "pose.hpp"

#include <ipc/geometry/normal.hpp>
#include <ipc/math/sinc.hpp>

namespace ipc::rigid {

Eigen::Matrix3d
rotation_vector_to_matrix(Eigen::ConstRef<Eigen::Vector3d> theta)
{
    const double sinc_angle = sinc_norm_x<double>(theta);
    const double sinc_half_angle = sinc_norm_x<double>((theta / 2).eval());
    const Eigen::Matrix3d K = cross_product_matrix(theta);
    const Eigen::Matrix3d K2 = K * K;
    Eigen::Matrix3d R =
        sinc_angle * K + 0.5 * sinc_half_angle * sinc_half_angle * K2;
    R.diagonal().array() += 1.0;
    return R;
}

Eigen::Matrix<double, 9, 3>
rotation_vector_to_matrix_jacobian(Eigen::ConstRef<Eigen::Vector3d> theta)
{
    assert(theta.size() == 3); // Only valid for 3D rotation vectors

    const double s = sinc_norm_x<double>(theta);
    const double sh = sinc_norm_x<double>((theta / 2).eval());
    const double sh_sq = sh * sh;
    const Eigen::Matrix3d K = cross_product_matrix(theta);

    const Eigen::Vector3d ds = sinc_norm_x_grad(theta);
    const Eigen::Vector3d dsh = 0.5 * sinc_norm_x_grad((theta / 2).eval());

    std::array<Eigen::Matrix3d, 3> dK;
    dK[0] << 0, 0, 0, 0, 0, -1, 0, 1, 0;
    dK[1] << 0, 0, 1, 0, 0, 0, -1, 0, 0;
    dK[2] << 0, -1, 0, 1, 0, 0, 0, 0, 0;

    std::array<Eigen::Matrix3d, 3> s_dK;
    s_dK[0] << 0, 0, 0, 0, 0, -s, 0, s, 0;
    s_dK[1] << 0, 0, s, 0, 0, 0, -s, 0, 0;
    s_dK[2] << 0, -s, 0, s, 0, 0, 0, 0, 0;

    const Eigen::Matrix3d sh_K2 = sh * K * K;

    Eigen::Matrix<double, 9, 3> dR;
    for (int i = 0; i < 3; ++i) {
        const Eigen::Matrix3d K_dKi = K * dK[i];
        dR.col(i) = (ds(i) * K + s_dK[i] + dsh(i) * sh_K2
                     + (0.5 * sh_sq) * (K_dKi + K_dKi.transpose()))
                        .reshaped();
    }

    return dR;
}

Eigen::Matrix<double, 9, 9>
rotation_vector_to_matrix_hessian(Eigen::ConstRef<Eigen::Vector3d> theta)
{
    assert(theta.size() == 3); // Only valid for 3D rotation vectors

    const double sh = sinc_norm_x<double>(theta / 2);
    const double sh_sq = sh * sh;
    const Eigen::Matrix3d K = cross_product_matrix(theta);

    const Eigen::Vector3d ds = sinc_norm_x_grad(theta);
    const Eigen::Vector3d dsh = 0.5 * sinc_norm_x_grad(theta / 2);
    const Eigen::Matrix3d d2s = sinc_norm_x_hess(theta);
    const Eigen::Matrix3d d2sh = 0.25 * sinc_norm_x_hess(theta / 2);

    std::array<Eigen::Matrix3d, 3> dK;
    dK[0] << 0, 0, 0, 0, 0, -1, 0, 1, 0;
    dK[1] << 0, 0, 1, 0, 0, 0, -1, 0, 0;
    dK[2] << 0, -1, 0, 1, 0, 0, 0, 0, 0;

    // NOTE: ∂²K/∂θ² = 0

    const Eigen::Matrix3d K2 = K * K;

    std::array<Eigen::Matrix3d, 3> K_dK;
    K_dK[0] = K * dK[0];
    K_dK[1] = K * dK[1];
    K_dK[2] = K * dK[2];

    std::array<Eigen::Matrix3d, 3> K_dK_plus_K_dK_T;
    K_dK_plus_K_dK_T[0] = K_dK[0] + K_dK[0].transpose();
    K_dK_plus_K_dK_T[1] = K_dK[1] + K_dK[1].transpose();
    K_dK_plus_K_dK_T[2] = K_dK[2] + K_dK[2].transpose();

    // Flatten in column-major order
    Eigen::Matrix<double, 9, 9> d2R;
    for (int j = 0; j < 3; ++j) {
        const Eigen::Matrix3d K_dKj = K * dK[j];
        for (int i = j; i < 3; ++i) {
            const Eigen::Matrix3d dKj_dKi = dK[j] * dK[i];
            d2R.col(j * 3 + i) =
                (d2s(i, j) * K + ds(i) * dK[j] + ds(j) * dK[i]
                 + (d2sh(i, j) * sh + dsh(i) * dsh(j)) * K2
                 + (sh * dsh(i)) * K_dK_plus_K_dK_T[j]
                 + (sh * dsh(j)) * K_dK_plus_K_dK_T[i]
                 + (0.5 * sh_sq) * (dKj_dKi + dKj_dKi.transpose()))
                    .reshaped();
            if (i != j) {
                d2R.col(i * 3 + j) = d2R.col(j * 3 + i); // Symmetric
            }
        }
    }
    return d2R;
}

Eigen::Vector3d rotation_matrix_to_vector(Eigen::ConstRef<Eigen::Matrix3d> R)
{
#if false // NOLINT
    // Eigen does this conversion by going from SO(3) -> Quaternion -> 𝔰𝔬(3),
    // but we can do it directly from SO(3) -> 𝔰𝔬(3). In random benchmarking,
    // this is about 2x faster than the Eigen implementation. However, this
    // might be less accurate for some edge cases as it requires a sin and acos
    // where as Eigen's approach uses atan2 and sqrt.

    assert(R.isUnitary(1e-9)); // Ensure it's a rotation matrix

    Eigen::Vector3d r;
    // θ ∈ [0, π)
    double theta = acos((std::clamp(R.trace(), -1.0, 3.0) - 1) / 2.0);
    assert(std::isfinite(theta));

    // R = I + sin(θ)K + (1 - cos(θ))K² where K is the cross-product matrix of
    // r/θ
    if (theta <= 1e-5) {
        // θ ≈ 0 ⟹ R ≈ I + θK
        // No need to divide by sin(θ)≈θ since the off-diagonals are ±rᵢθ
        r.x() = (R(2, 1) - R(1, 2)) / 2;
        r.y() = (R(0, 2) - R(2, 0)) / 2;
        r.z() = (R(1, 0) - R(0, 1)) / 2;
        // Θ is already part of the off-diagonal elements, so no need to scale
    } else if (theta >= M_PI - 1e-5) {
        // θ ≈ π ⟹ R ≈ I - (θ - π)K + 2K²
        // Ensure positive for sqrt
        assert(R(0, 0) > -1);
        r.x() = sqrt((R(0, 0) + 1) / 2);
        // Determine the other components based on the first component
        // R₀,₁ = 2r₀r₁+(θ-π)r₂, R₀,₂ = 2r₀r₂-(θ-π)r₁, R₁,₂ = 2r₁r₂+(θ-π)r₀
        r.y() = R(0, 1) / (2 * r.x());
        r.z() = R(0, 2) / (2 * r.x());
        // r.normalize();      // Normalize just in case
        r.array() *= theta; // Scale by the angle
    } else {
        const double sin_theta = sin(theta);
        r.x() = (R(2, 1) - R(1, 2)) / (2 * sin_theta);
        r.y() = (R(0, 2) - R(2, 0)) / (2 * sin_theta);
        r.z() = (R(1, 0) - R(0, 1)) / (2 * sin_theta);
        // r.normalize(); // Normalize just in case
        r *= theta; // Scale by the angle
    }
    assert(r.array().isFinite().all());

#ifndef NDEBUG
    Eigen::AngleAxisd eigen(R);
    assert(
        (eigen.angle() * eigen.axis()).isApprox(r, 1e-4)
        || (theta >= M_PI - 1e-5 // Eigen can flip the sign of angles close to π
            && (-eigen.angle() * eigen.axis()).isApprox(r, 1e-4)));
#endif

    return r;
#else
    // Use Eigen's built-in function for conversion
    Eigen::AngleAxisd angle_axis(R);
    return angle_axis.angle() * angle_axis.axis();
#endif
}

Eigen::MatrixXd
Pose::transform_vertices_jacobian(Eigen::ConstRef<Eigen::MatrixXd> V) const
{
    const int dim = V.cols();
    const int ndof = position.size() + rotation.size();
    assert(dim == position.size());

    Eigen::MatrixXd jac = Eigen::MatrixXd::Zero(V.size(), ndof);

    // Derivative w.r.t. position is identity
    for (int d = 0; d < dim; ++d) {
        jac.block(d * V.rows(), d, V.rows(), 1).setOnes();
    }

    // Precompute dR/dθ
    MatrixMax<double, 9, 3> dR_dtheta;
    if (dim == 2) {
        dR_dtheta.resize(4, 1);
        // 2D rotation matrix:
        // R(θ) = [cos(θ) -sin(θ);
        //         sin(θ)  cos(θ)]
        const double dcos = -std::sin(rotation(0));
        const double dsin = std::cos(rotation(0));
        dR_dtheta(0, 0) = dcos;
        dR_dtheta(1, 0) = dsin;
        dR_dtheta(2, 0) = -dsin;
        dR_dtheta(3, 0) = dcos;
    } else {
        dR_dtheta = rotation_vector_to_matrix_jacobian(rotation);
    }

    // Derivative w.r.t. rotation
    for (int j = 0; j < rotation.size(); ++j) {
        jac.block(0, dim + j, V.size(), 1) =
            (V * dR_dtheta.col(j).reshaped(dim, dim).transpose()).reshaped();
    }

    return jac;
}

Eigen::MatrixXd
Pose::transform_vertices_hessian(Eigen::ConstRef<Eigen::MatrixXd> V) const
{
    const int dim = V.cols();
    const int ndof = position.size() + rotation.size();
    assert(dim == position.size());

    Eigen::MatrixXd hess = Eigen::MatrixXd::Zero(V.size(), ndof * ndof);

    // Derivative of position Jacobian is zero, so only need to compute
    // derivatives involving rotation.

    // Precompute d^2R/dθ^2
    MatrixMax<double, 9, 9> d2R_dtheta2;
    if (dim == 2) {
        d2R_dtheta2.resize(4, 1);
        // 2D rotation matrix:
        // R(θ) = [cos(θ) -sin(θ);
        //         sin(θ)  cos(θ)]
        const double d2cos = -std::cos(rotation(0));
        const double d2sin = -std::sin(rotation(0));
        d2R_dtheta2(0, 0) = d2cos;
        d2R_dtheta2(1, 0) = d2sin;
        d2R_dtheta2(2, 0) = -d2sin;
        d2R_dtheta2(3, 0) = d2cos;
    } else {
        d2R_dtheta2 = rotation_vector_to_matrix_hessian(rotation);
    }

    // Derivative w.r.t. rotation
    for (int j = position.size(); j < ndof; ++j) {
        for (int i = position.size(); i < ndof; ++i) {
            const int k = (j - dim) * rotation.size() + (i - dim);
            hess.col(j * ndof + i) =
                (V * d2R_dtheta2.col(k).reshaped(dim, dim).transpose())
                    .reshaped();
        }
    }

    return hess;
}

} // namespace ipc::rigid