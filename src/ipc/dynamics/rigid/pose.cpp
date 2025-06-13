#include "pose.hpp"

#include <ipc/utils/sinc.hpp>

namespace ipc::rigid {

namespace {
    inline Eigen::Matrix3d
    cross_product_matrix(Eigen::ConstRef<Eigen::Vector3d> x)
    {
        Eigen::Matrix3d X;
        X << 0, -x.z(), x.y(), //
            x.z(), 0, -x.x(),  //
            -x.y(), x.x(), 0;
        return X;
    }
} // namespace

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

Eigen::Matrix<double, 3, 9>
rotation_vector_to_matrix_jacobian(Eigen::ConstRef<Eigen::Vector3d> theta)
{
    assert(theta.size() == 3); // Only valid for 3D rotation vectors

    const double s = sinc_norm_x<double>(theta);
    const double sh = sinc_norm_x<double>(theta / 2);
    const double sh_sq = sh * sh;
    const Eigen::Matrix3d K = cross_product_matrix(theta);

    const Eigen::Vector3d ds = sinc_norm_x_grad(theta);
    const Eigen::Vector3d dsh = sinc_norm_x_grad(theta / 2);

    std::array<Eigen::Matrix3d, 3> dK;
    dK[0] << 0, 0, 0, 0, 0, -1, 0, 1, 0;
    dK[1] << 0, 0, 1, 0, 0, 0, -1, 0, 0;
    dK[2] << 0, -1, 0, 1, 0, 0, 0, 0, 0;

    Eigen::Matrix3d s_plus_sh_sq_K = sh_sq * K;
    s_plus_sh_sq_K.diagonal().array() += s;

    const Eigen::Matrix3d sh_K2 = sh * K * K;

    Eigen::Matrix<double, 3, 9> dR;
    for (int i = 0; i < 3; ++i) {
        dR.row(i) =
            (ds(i) * K + s_plus_sh_sq_K * dK[i] + dsh(i) * sh_K2).reshaped();
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
    const Eigen::Vector3d dsh = sinc_norm_x_grad(theta / 2);
    const Eigen::Matrix3d d2s = sinc_norm_x_hess(theta);
    const Eigen::Matrix3d d2sh = sinc_norm_x_hess(theta / 2);

    std::array<Eigen::Matrix3d, 3> dK;
    dK[0] << 0, 0, 0, 0, 0, -1, 0, 1, 0;
    dK[1] << 0, 0, 1, 0, 0, 0, -1, 0, 0;
    dK[2] << 0, -1, 0, 1, 0, 0, 0, 0, 0;

    // NOTE: ∂²K/∂θ² = 0

    const Eigen::Matrix3d sh_K2 = sh * K * K;

    Eigen::Matrix<double, 9, 9> d2R;
    for (int i = 0; i < 3; ++i) {
        const Eigen::Matrix3d two_sh_K_dKi = 2 * sh * K * dK[i];

        for (int j = 0; j < 3; ++j) {
            d2R.row(i * 3 + j) =
                (d2s(i, j) * K + ds(i) * dK[j] + ds(j) * dK[i]
                 + dsh(j) * two_sh_K_dKi + sh_sq * dK[j] * dK[i]
                 + d2sh(i, j) * sh_K2 + 2 * dsh(i) * K * dK[j])
                    .reshaped();
        }
    }
    return d2R;
}

Eigen::Vector3d rotation_matrix_to_vector(Eigen::Matrix3d R)
{
    Eigen::Vector3d r;
    assert(R.trace() >= -1 && R.trace() <= 3);  // Ensure acos is valid
    double theta = acos((R.trace() - 1) / 2.0); // θ ∈ [0, π)

    // R = I + sin(θ)K + (1 - cos(θ))K² where K is the cross-product matrix of
    // r/θ
    if (theta < 1e-6) {
        // θ ≈ 0 ⟹ R ≈ I + θK
        // No need to divide by sin(θ)≈θ since the off-diagonals are ±rᵢθ
        r.x() = (R(2, 1) - R(1, 2)) / 2;
        r.y() = (R(0, 2) - R(2, 0)) / 2;
        r.z() = (R(1, 0) - R(0, 1)) / 2;
        // Θ is already part of the off-diagonal elements, so no need to scale
    } else if (theta > M_PI - 1e-6) {
        // θ ≈ π ⟹ R ≈ I - (θ - π)K + 2K²
        // Ensure positive for sqrt
        assert(R(0, 0) > -1);
        r.x() = sqrt((R(0, 0) + 1) / 2);
        // Determine the other components based on the first component
        // R₀,₁ = 2r₀r₁+(θ-π)r₂, R₀,₂ = 2r₀r₂-(θ-π)r₁, R₁,₂ = 2r₁r₂+(θ-π)r₀
        r.y() = R(0, 1) / (2 * r.x());
        r.z() = R(0, 2) / (2 * r.x());
        r.normalize();      // Normalize just in case
        r.array() *= theta; // Scale by the angle
    } else {
        r.x() = (R(2, 1) - R(1, 2)) / (2 * sin(theta));
        r.y() = (R(0, 2) - R(2, 0)) / (2 * sin(theta));
        r.z() = (R(1, 0) - R(0, 1)) / (2 * sin(theta));
        r.normalize(); // Normalize just in case
        r *= theta;    // Scale by the angle
    }

    return r;
}

template struct Pose<RotationVector>;
template struct Pose<RotationMatrix>;

} // namespace ipc::rigid