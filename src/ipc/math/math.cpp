#include "math.hpp"

#include <ipc/gcp/common.hpp>
#include <ipc/geometry/normal.hpp>
#include <ipc/utils/autodiff_types.hpp>

#include <Eigen/SVD>

#include <algorithm>

namespace ipc {
namespace {

    // Rotation-variant SVD of a fixed-size rows x cols matrix A = U diag(sigma)
    // Vᵀ, where U is rows x rows, V is cols x cols, and sigma has
    // min(rows, cols) entries. We compute the full U and V, then adjust signs
    // so that both are proper rotations (det = +1) whenever we can.
    //
    // For a square matrix, the only freedom we have is the sign of the
    // smallest singular value. So we fold any reflection into it, following the
    // standard signed SVD. This means sigma can go negative when det(A) < 0.
    //
    // For a rectangular matrix, the larger factor has extra columns that span
    // the null space of A, and flipping one of them changes nothing about
    // U diag(sigma) Vᵀ. We use that spare column to fix the determinant of the
    // larger factor. For the smaller factor, we flip the same column in both
    // U and V, which also leaves the product untouched. The tradeoff is that
    // sigma stays non-negative here, so U_thin Vᵀ (the nearest matrix with
    // orthonormal columns or rows) is exactly the polar factor.
    template <int rows, int cols>
    void svd_rv(
        const Eigen::Matrix<double, rows, cols>& A,
        Eigen::Matrix<double, rows, rows>& U,
        Eigen::Matrix<double, cols, cols>& V,
        Eigen::Vector<double, std::min(rows, cols)>& sigma)
    {
        constexpr int k = std::min(rows, cols);
        // Compute-time U/V options (runtime options are deprecated in Eigen).
        constexpr int COMPUTE_FULL_UV =
            Eigen::ComputeFullU | Eigen::ComputeFullV;

        Eigen::JacobiSVD<Eigen::Matrix<double, rows, cols>, COMPUTE_FULL_UV>
            svd(A);

        U = svd.matrixU();
        V = svd.matrixV();
        sigma = svd.singularValues();

        if constexpr (rows == cols) {
            using MatrixType = Eigen::Matrix<double, rows, rows>;
            MatrixType L = MatrixType::Identity();
            L(rows - 1, rows - 1) = (U * V.transpose()).determinant();

            const double detU = U.determinant();
            const double detV = V.determinant();
            if (detU < 0 && detV > 0) {
                U = U * L;
            } else if (detU > 0 && detV < 0) {
                V = V * L;
            }
            sigma = L * sigma;
        } else {
            // Flipping column k-1 of both U and V keeps U diag(sigma) Vᵀ the
            // same, and flipping a null-space column of the larger factor
            // changes nothing at all.
            const auto fix_signs = [](auto& small, auto& big) {
                if (small.determinant() < 0) {
                    small.col(k - 1) *= -1;
                    big.col(k - 1) *= -1;
                }
                if (big.determinant() < 0) {
                    big.col(big.cols() - 1) *= -1;
                }
            };
            if constexpr (rows > cols) {
                fix_signs(V, U);
            } else {
                fix_signs(U, V);
            }
        }
    }

    // The nearest rotation (rows x cols) built from the rotation-variant SVD.
    template <int rows, int cols>
    MatrixMax3d nearest_rotation_fixed(Eigen::ConstRef<MatrixMax3d> A)
    {
        constexpr int k = std::min(rows, cols);
        Eigen::Matrix<double, rows, rows> U;
        Eigen::Matrix<double, cols, cols> V;
        Eigen::Vector<double, k> sigma;
        svd_rv<rows, cols>(Eigen::Matrix<double, rows, cols>(A), U, V, sigma);
        return U.template leftCols<k>() * V.template leftCols<k>().transpose();
    }

    /*
        Mathematica script

        t1x = Table[t1[k], {k, 0, 2}];
        t2x = Table[t2[k], {k, 0, 2}];
        t = Flatten[{t1x, t2x}];
        n = Cross[t1x, t2x];
        CForm[Sum[ddy[i - 1, j - 1] Transpose[{D[n[[j]], {t}]}].{D[n[[i]], {t}]}
       + dy[i - 1] D[n[[i]], {t}, {t}], {i, 1, 3}, {j, 1, 3}]]
    */

} // namespace

MatrixMax3d nearest_rotation(Eigen::ConstRef<MatrixMax3d> A)
{
    assert(A.rows() >= 2 && A.rows() <= 3 && A.cols() >= 2 && A.cols() <= 3);
    if (A.rows() == 2 && A.cols() == 2) {
        return nearest_rotation_fixed<2, 2>(A);
    } else if (A.rows() == 3 && A.cols() == 3) {
        return nearest_rotation_fixed<3, 3>(A);
    } else if (A.rows() == 3 && A.cols() == 2) {
        return nearest_rotation_fixed<3, 2>(A);
    } else {
        return nearest_rotation_fixed<2, 3>(A);
    }
}

HeavisideType OrientationTypes::compute_type(
    const double val, const double alpha, const double beta)
{
    if (val <= -alpha) {
        return HeavisideType::ZERO;
    }
    if (val >= beta) {
        return HeavisideType::ONE;
    }

    return HeavisideType::VARIANT;
}

void OrientationTypes::set_size(const int size)
{
    m_size = size;
    tangent_types.assign(m_size, HeavisideType::VARIANT);
    normal_types.assign(m_size, HeavisideType::VARIANT);
}

double opposite_direction_penalty(
    Eigen::ConstRef<Eigen::Vector3d> t,
    Eigen::ConstRef<Eigen::Vector3d> d,
    const double alpha,
    const double beta)
{
    return Math<double>::smooth_heaviside(d.dot(t) / t.norm(), alpha, beta);
}

GradientType<6> opposite_direction_penalty_grad(
    Eigen::ConstRef<Eigen::Vector3d> t,
    Eigen::ConstRef<Eigen::Vector3d> d,
    const double alpha,
    const double beta)
{
    auto [tn, tn_grad] = normalization_and_jacobian(t);
    const double a = d.dot(tn);
    const double y = Math<double>::smooth_heaviside(a, alpha, beta);
    const double dy = Math<double>::smooth_heaviside_grad(a, alpha, beta);

    Vector6d grad;
    grad << tn_grad * dy * d, dy * tn;
    return std::make_tuple(y, grad);
}

HessianType<6> opposite_direction_penalty_hess(
    Eigen::ConstRef<Eigen::Vector3d> t,
    Eigen::ConstRef<Eigen::Vector3d> d,
    const double alpha,
    const double beta)
{
    auto [tn, tn_grad, tn_hess] = normalization_and_jacobian_and_hessian(t);
    const double a = d.dot(tn);
    const double y = Math<double>::smooth_heaviside(a, alpha, beta);
    const double dy = Math<double>::smooth_heaviside_grad(a, alpha, beta);
    const double ddy = Math<double>::smooth_heaviside_hess(a, alpha, beta);

    Vector6d grad;
    grad << dy * d, dy * tn;

    Matrix6d hess;
    hess << d * ddy * d.transpose(),
        d * ddy * tn.transpose() + Eigen::Matrix3d::Identity() * dy,
        tn * ddy * d.transpose() + Eigen::Matrix3d::Identity() * dy,
        tn * ddy * tn.transpose();

    // chain rule of vector t normalize
    hess.topRows(3) = tn_grad * hess.topRows(3);
    hess.leftCols(3) = hess.leftCols(3) * tn_grad;
    hess.topLeftCorner(3, 3) +=
        grad(0) * tn_hess[0] + grad(1) * tn_hess[1] + grad(2) * tn_hess[2];

    grad.head(3) = tn_grad * grad.head(3);

    return std::make_tuple(y, grad, hess);
}

// assume unit vector d
double negative_orientation_penalty(
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2,
    Eigen::ConstRef<Eigen::Vector3d> d,
    const double alpha,
    const double beta)
{
    const Eigen::Vector3d n = t1.cross(t2);
    return opposite_direction_penalty(n, d, alpha, beta);
}

GradientType<9> negative_orientation_penalty_grad(
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2,
    Eigen::ConstRef<Eigen::Vector3d> d,
    const double alpha,
    const double beta)
{
    const Eigen::Vector3d n = t1.cross(t2);
    auto [y, dy] = opposite_direction_penalty_grad(n, d, alpha, beta);

    /*
        Mathematica script

        t1 = {t10, t11, t12};
        t2 = {t20, t21, t22};
        n = Cross[t1, t2];
        Sum[dy[k - 1] D[n[[k]], {{t1, t2}}], {k, 1, 3}] // MatrixForm
    */

    Vector9d grad;
    grad << -t2(2) * dy(1) + t2(1) * dy(2), -t2(0) * dy(2) + t2(2) * dy(0),
        -t2(1) * dy(0) + t2(0) * dy(1), -t1(1) * dy(2) + t1(2) * dy(1),
        -t1(2) * dy(0) + t1(0) * dy(2), -t1(0) * dy(1) + t1(1) * dy(0),
        dy.tail(3);
    return std::make_tuple(y, grad);
}

HessianType<9> negative_orientation_penalty_hess(
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2,
    Eigen::ConstRef<Eigen::Vector3d> d,
    const double alpha,
    const double beta)
{
    const Eigen::Vector3d n = t1.cross(t2);
    const Eigen::Matrix<double, 3, 6> cross_grad =
        cross_product_gradient(t1, t2);
    const std::array<Matrix6d, 3> cross_hess = cross_product_hessian(t1, t2);
    auto [y, dy, ddy] = opposite_direction_penalty_hess(n, d, alpha, beta);

    Vector9d grad;
    grad.tail(3) = dy.tail(3);
    grad.head(6) = cross_grad.transpose() * dy.head(3);

    Matrix9d hess;
    hess.bottomRightCorner(3, 3) = ddy.bottomRightCorner(3, 3);
    hess.topLeftCorner(6, 6) =
        cross_grad.transpose() * ddy.topLeftCorner(3, 3) * cross_grad
        + cross_hess[0] * dy(0) + cross_hess[1] * dy(1) + cross_hess[2] * dy(2);
    hess.topRightCorner(6, 3) =
        cross_grad.transpose() * ddy.topRightCorner(3, 3);
    hess.bottomLeftCorner(3, 6) = ddy.bottomLeftCorner(3, 3) * cross_grad;

    return std::make_tuple(y, grad, hess);
}

Eigen::Matrix<double, 3, 6> cross_product_gradient(
    Eigen::ConstRef<Eigen::Vector3d> t1, Eigen::ConstRef<Eigen::Vector3d> t2)
{
    Eigen::Matrix<double, 3, 6> grad;
    grad << 0, t2(2), -t2(1), 0, -t1(2), t1(1), -t2(2), 0, t2(0), t1(2), 0,
        -t1(0), t2(1), -t2(0), 0, -t1(1), t1(0), 0;

    return grad;
}

std::array<Matrix6d, 3> cross_product_hessian(
    Eigen::ConstRef<Eigen::Vector3d> t1, Eigen::ConstRef<Eigen::Vector3d> t2)
{
    std::array<Matrix6d, 3> hess;
    hess.fill(Matrix6d::Zero());
    hess[0](1, 5) = 1;
    hess[0](5, 1) = 1;
    hess[0](2, 4) = -1;
    hess[0](4, 2) = -1;

    hess[1](0, 5) = -1;
    hess[1](5, 0) = -1;
    hess[1](2, 3) = 1;
    hess[1](3, 2) = 1;

    hess[2](0, 4) = 1;
    hess[2](4, 0) = 1;
    hess[2](1, 3) = -1;
    hess[2](3, 1) = -1;

    return hess;
}

} // namespace ipc