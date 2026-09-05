#include "angle.hpp"

#include <ipc/config.hpp>
#include <ipc/geometry/normal.hpp>
#include <ipc/math/scalar_math.hpp>
#include <ipc/utils/simd.hpp>

#include <array>
#include <utility>

namespace ipc::detail {

namespace {
    /// @brief Jacobians of the two triangle normals with respect to all 12 DOFs
    /// (x0, x1, x2, x3).
    ///
    /// n0 = normal(x0, x1, x2) does not depend on x3, so its last block is
    /// zero. n1 = normal(x1, x0, x3) comes back in (x1, x0, x3) order, so we
    /// permute its columns into place and zero the x2 block.
    template <typename T>
    inline std::pair<Eigen::Matrix<T, 3, 12>, Eigen::Matrix<T, 3, 12>>
    dihedral_normal_jacobians(
        Eigen::ConstRef<Eigen::Vector3<T>> x0,
        Eigen::ConstRef<Eigen::Vector3<T>> x1,
        Eigen::ConstRef<Eigen::Vector3<T>> x2,
        Eigen::ConstRef<Eigen::Vector3<T>> x3)
    {
        Eigen::Matrix<T, 3, 12> dn0_dx;
        dn0_dx.template leftCols<9>() = triangle_normal_jacobian(x0, x1, x2);
        dn0_dx.template rightCols<3>().setZero();

        Eigen::Matrix<T, 3, 12> dn1_dx;
        const std::array<int, 9> idx = { { 3, 4, 5, 0, 1, 2, 9, 10, 11 } };
        dn1_dx(Eigen::all, idx) = triangle_normal_jacobian(x1, x0, x3);
        dn1_dx.template middleCols<3>(6).setZero();

        return { dn0_dx, dn1_dx };
    }
} // namespace

template <typename T>
T dihedral_angle(
    Eigen::ConstRef<Eigen::Vector3<T>> x0,
    Eigen::ConstRef<Eigen::Vector3<T>> x1,
    Eigen::ConstRef<Eigen::Vector3<T>> x2,
    Eigen::ConstRef<Eigen::Vector3<T>> x3)
{
    const Eigen::Vector3<T> n0 = triangle_normal(x0, x1, x2);
    const Eigen::Vector3<T> n1 = triangle_normal(x1, x0, x3);
    const Eigen::Vector3<T> e = normalized(x1 - x0);

    const T sin_theta = n0.cross(n1).dot(e);
    const T cos_theta = n0.dot(n1);

    return ipc::atan2(sin_theta, cos_theta);
}

template <typename T>
Eigen::Vector<T, 12> dihedral_angle_gradient(
    Eigen::ConstRef<Eigen::Vector3<T>> x0,
    Eigen::ConstRef<Eigen::Vector3<T>> x1,
    Eigen::ConstRef<Eigen::Vector3<T>> x2,
    Eigen::ConstRef<Eigen::Vector3<T>> x3)
{
    const Eigen::Vector3<T> n0 = triangle_normal(x0, x1, x2);
    const Eigen::Vector3<T> n1 = triangle_normal(x1, x0, x3);
    const Eigen::Vector3<T> e = normalized(x1 - x0);

    // --- Normal gradients ---

    const auto [dn0_dx, dn1_dx] = dihedral_normal_jacobians<T>(x0, x1, x2, x3);

    // --- Angle gradient ---

    const Eigen::Vector<T, 12> dcos_dx =
        dn0_dx.transpose() * n1 + dn1_dx.transpose() * n0;

    const Eigen::Vector<T, 12> dsin_dx =
        (cross_product_matrix<T>(n0) * dn1_dx
         - cross_product_matrix<T>(n1) * dn0_dx)
            .transpose()
        * e;

    // --- Product rule ---

    const T sin_theta = n0.cross(n1).dot(e);
    const T cos_theta = n0.dot(n1);

    return dsin_dx * cos_theta - dcos_dx * sin_theta;
}

template <typename T>
Eigen::Matrix<T, 12, 12> dihedral_angle_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> x0,
    Eigen::ConstRef<Eigen::Vector3<T>> x1,
    Eigen::ConstRef<Eigen::Vector3<T>> x2,
    Eigen::ConstRef<Eigen::Vector3<T>> x3)
{
    using Matrix12 = Eigen::Matrix<T, 12, 12>;

    const Eigen::Vector3<T> n0 = triangle_normal(x0, x1, x2);
    const Eigen::Vector3<T> n1 = triangle_normal(x1, x0, x3);

    // -------------------------------------------------------------------------
    // Jacobian of n0 and n1 w.r.t. all 12 DOFs
    // -------------------------------------------------------------------------

    const auto [dn0_dx, dn1_dx] = dihedral_normal_jacobians<T>(x0, x1, x2, x3);

    // -------------------------------------------------------------------------
    // Hessians of n0 and n1 — shape (27 × 9) each, then embedded in 12-DOF
    // space
    //
    // dn0: lives in the first 9 DOFs (x0, x1, x2).
    //   d2n0_full[3*i, j]  with i,j ∈ {0..8}
    //
    // dn1 = triangle_normal(x1, x0, x3): lives in DOFs {x1, x0, x3}
    //   mapped back: local order (x1→3..5, x0→0..2, x3→9..11)
    // -------------------------------------------------------------------------

    // Raw hessians in local (9-DOF) coordinate systems
    const Eigen::Matrix<T, 27, 9> d2n0_local =
        triangle_normal_hessian(x0, x1, x2);
    const Eigen::Matrix<T, 27, 9> d2n1_local =
        triangle_normal_hessian(x1, x0, x3);

    // We will access them on-the-fly via the index maps rather than building
    // full (12-DOF) tensors.  The local DOF → global DOF maps are:
    //   n0: local col k → global col k  (for k in 0..8), col 9..11 → zero
    //   n1: local col k → global col idx[k]  (idx = {3,4,5,0,1,2,9,10,11})

    // Helper: retrieve d2nX_k / (d xi d xj) for normal component k (0..2),
    //         global DOF indices p and q.
    // Returns 0 if the global index falls outside the local support.

    // Global to local maps  (-1 means "not in support → zero")
    // n0 support: global 0..8  → local 0..8;  global 9..11 → -1
    auto g2l_n0 = [](int g) -> int { return (g < 9) ? g : -1; };

    // n1 support: global {0,1,2}→local{3,4,5},  global{3,4,5}→local{0,1,2},
    //             global{9,10,11}→local{6,7,8},  global{6,7,8}→-1
    const std::array<int, 12> g2l_n1_map = { { 3, 4, 5, 0, 1, 2, -1, -1, -1, 6,
                                               7, 8 } };
    auto g2l_n1 = [&g2l_n1_map](int g) -> int { return g2l_n1_map[g]; };

    // -------------------------------------------------------------------------
    // Hessian of the normalized edge e — non-zero only for global DOFs 0..5.
    // We store two 3×3 blocks:
    //   He_block[k](p,q) = ∂²e_k / (∂x_{p} ∂x_{q})   with x_p, x_q ∈ {x0,x1}
    //
    // e = f(v) = v/||v||,  v = x1 - x0
    // ∂e_k/∂v_p = J_e[k,p]
    // ∂²e_k/(∂v_p ∂v_q) = He[k][p,q]   (from
    // normalization_and_jacobian_and_hessian)
    //
    // Chain rule for v = x1 - x0:
    //   ∂v/∂x0 = -I,  ∂v/∂x1 = +I
    // ∂²e_k/(∂x0_p ∂x0_q) = +He[k][p,q]
    // ∂²e_k/(∂x1_p ∂x1_q) = +He[k][p,q]
    // ∂²e_k/(∂x0_p ∂x1_q) = -He[k][p,q]
    // ∂²e_k/(∂x1_p ∂x0_q) = -He[k][p,q]
    // -------------------------------------------------------------------------

    const auto [e, Je, He] =
        normalization_and_jacobian_and_hessian<T, 3>(x1 - x0);

    // He is std::array<Eigen::Matrix3<T>, 3> where He[k] is the 3×3 Hessian for
    // component k of the normalized vector.

    // -------------------------------------------------------------------------
    // Jacobian of the normalized edge direction e = (x1 - x0) / ||x1 - x0||
    // Only x0 (cols 0..2) and x1 (cols 3..5) contribute.
    // ∂e/∂x0 = -J_norm,  ∂e/∂x1 = +J_norm
    // -------------------------------------------------------------------------

    Eigen::Matrix<T, 3, 12> de_dx = Eigen::Matrix<T, 3, 12>::Zero();
    de_dx.template leftCols<3>() = -Je;   // ∂e/∂x0
    de_dx.template middleCols<3>(3) = Je; // ∂e/∂x1

    // -------------------------------------------------------------------------
    // Scalars  s = (n0 × n1) · e  and  c = n0 · n1
    // -------------------------------------------------------------------------

    const Eigen::Vector3<T> m = n0.cross(n1); // n0 × n1
    const T sin_theta = m.dot(e);
    const T cos_theta = n0.dot(n1);

    // -------------------------------------------------------------------------
    // Jacobian of s and c (12-vectors) — re-derived consistently with gradient
    // -------------------------------------------------------------------------

    // dm_dx = ∂(n0×n1)/∂x  (3×12)
    const Eigen::Matrix<T, 3, 12> dm_dx = cross_product_matrix<T>(n0) * dn1_dx
        - cross_product_matrix<T>(n1) * dn0_dx;

    const Eigen::Vector<T, 12> dcos_dx =
        dn0_dx.transpose() * n1 + dn1_dx.transpose() * n0;

    const Eigen::Vector<T, 12> dsin_dx =
        dm_dx.transpose() * e + de_dx.transpose() * m;

    const Eigen::Vector<T, 12> dtheta_dx =
        dsin_dx * cos_theta - dcos_dx * sin_theta;

    // -------------------------------------------------------------------------
    // Hessian of c = n0 · n1
    //
    // ∂²c/(∂xp ∂xq) = (∂²n0/(∂xp ∂xq)) · n1
    //                + (∂n0/∂xp) · (∂n1/∂xq)
    //                + (∂n1/∂xp) · (∂n0/∂xq)
    //                + n0 · (∂²n1/(∂xp ∂xq))
    // -------------------------------------------------------------------------

    Matrix12 H_cos;

    // Cross-Jacobian terms (symmetric)
    H_cos = dn0_dx.transpose() * dn1_dx;
    H_cos += dn1_dx.transpose() * dn0_dx;

    // Second-derivative terms: sum over normal components k
    for (int k = 0; k < 3; k++) {
        // Contribution from n0's hessian contracted with n1[k]
        // and from n1's hessian contracted with n0[k]
        for (int p = 0; p < 12; p++) {
            const int lp0 = g2l_n0(p);
            const int lp1 = g2l_n1(p);
            for (int q = 0; q < 12; q++) {
                const int lq0 = g2l_n0(q);
                const int lq1 = g2l_n1(q);
                // d2n0_k / (dxp dxq)
                if (lp0 >= 0 && lq0 >= 0) {
                    H_cos(p, q) += n1(k) * d2n0_local(3 * lp0 + k, lq0);
                }
                // d2n1_k / (dxp dxq)
                if (lp1 >= 0 && lq1 >= 0) {
                    H_cos(p, q) += n0(k) * d2n1_local(3 * lp1 + k, lq1);
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Hessian of s = m · e = (n0 × n1) · e
    //
    // ∂²s/(∂xp ∂xq) = [∂²m/(∂xp ∂xq)] · e
    //                + [∂m/∂xp] · [∂e/∂xq]
    //                + [∂e/∂xp] · [∂m/∂xq]
    //                + m · [∂²e/(∂xp ∂xq)]
    // -------------------------------------------------------------------------

    Matrix12 H_sin;

    // Cross-Jacobian terms (m–e interaction)
    H_sin = dm_dx.transpose() * de_dx;
    H_sin += de_dx.transpose() * dm_dx;

    // Second derivative of e contracted with m — only DOFs 0..5 contribute
    // ∂²e_k/(∂v_p ∂v_q) with signs for x0/x1 blocks
    // Block (x0,x0): +He[k], Block (x1,x1): +He[k],
    // Block (x0,x1): -He[k], Block (x1,x0): -He[k]
    for (int k = 0; k < 3; k++) {
        const Eigen::Matrix3<T>& Hek = He[k]; // 3×3
        const T mk = m(k);
        // (x0, x0)
        H_sin.template block<3, 3>(0, 0) += mk * Hek;
        // (x1, x1)
        H_sin.template block<3, 3>(3, 3) += mk * Hek;
        // (x0, x1)
        H_sin.template block<3, 3>(0, 3) -= mk * Hek;
        // (x1, x0)
        H_sin.template block<3, 3>(3, 0) -= mk * Hek;
    }

    // Second derivative of m = n0 × n1 contracted with e:
    //
    // ∂²(n0 × n1)/(∂xp ∂xq)
    //   = -[∂n1/∂xq ×] · ∂n0/∂xp  - [n1×] · ∂²n0/(∂xp ∂xq)
    //   + [∂n0/∂xq ×] · ∂n1/∂xp   + [n0×] · ∂²n1/(∂xp ∂xq)
    //
    // Contracted with e:
    //   e · ∂²m/(∂xp ∂xq)
    //   = -eᵀ [∂n1/∂xq ×] ∂n0/∂xp  - eᵀ [n1×] ∂²n0/(∂xp ∂xq)
    //   + eᵀ [∂n0/∂xq ×] ∂n1/∂xp   + eᵀ [n0×] ∂²n1/(∂xp ∂xq)
    //
    // Note: eᵀ [v×] w = e · (v × w) = (e × v) · w  →  (-v × e) · w
    //   so  eᵀ [v×] = (-v × e)ᵀ = (e × v)ᵀ  (as a row vector)

    // Pre-compute (e × n0) for efficiency
    const Eigen::Vector3<T> e_cross_n0 = e.cross(n0);

    // --- Terms involving first × first derivatives (cross of Jacobian cols)
    // --- -eᵀ [dn1/dxq ×] dn0/dxp = (e × dn1/dxq) · dn0/dxp
    //                          = (e_cross_dn1_q) · dn0_dxp
    // +eᵀ [dn0/dxq ×] dn1/dxp = -(e × dn0/dxq) · dn1/dxp
    //    Wait: eᵀ [v×] w = e · (v × w). Let v = dn0/dxq, w = dn1/dxp
    //    => e · (dn0/dxq × dn1/dxp)
    //    Also: e · (v × w) = -(e × v) · w... let's just use dot directly.
    for (int q = 0; q < 12; q++) {
        const Eigen::Vector3<T> dn1_q = dn1_dx.col(q);
        const Eigen::Vector3<T> dn0_q = dn0_dx.col(q);
        // e · (dn1_q × dn0_p) for all p  →  negate the cross and dot with e
        // Term: -eᵀ [dn1_q×] dn0_dxp = e · (dn1_q × dn0_dxp) ... but
        //   [v×]w = v×w, so eᵀ [v×] w = e·(v×w) = (e×v)·w
        // -eᵀ [dn1_q ×] dn0_dxp = -(e × dn1_q) · dn0_dxp
        const Eigen::Vector3<T> neg_e_cross_dn1_q = -(e.cross(dn1_q));
        // +eᵀ [dn0_q ×] dn1_dxp = (e × dn0_q) · dn1_dxp
        const Eigen::Vector3<T> e_cross_dn0_q = e.cross(dn0_q);
        for (int p = 0; p < 12; p++) {
            H_sin(p, q) += neg_e_cross_dn1_q.dot(dn0_dx.col(p));
            H_sin(p, q) += e_cross_dn0_q.dot(dn1_dx.col(p));
        }
    }

    // --- Terms involving second derivatives of n0 and n1 contracted with e ---
    // -eᵀ [n1×] ∂²n0/(∂xp ∂xq)  = -(e × n1) · ∂²n0/(∂xp ∂xq)
    //                             = e_cross_n1_neg · ∂²n0/(∂xp ∂xq)
    // +eᵀ [n0×] ∂²n1/(∂xp ∂xq)  = (e × n0)      · ∂²n1/(∂xp ∂xq)
    //   where eᵀ [n1×] w = (e×n1)·w
    //         eᵀ [n0×] w = (e×n0)·w
    const Eigen::Vector3<T> neg_e_cross_n1 = -e.cross(n1);

    for (int p = 0; p < 12; p++) {
        const int lp0 = g2l_n0(p);
        const int lp1 = g2l_n1(p);
        for (int q = 0; q < 12; q++) {
            const int lq0 = g2l_n0(q);
            const int lq1 = g2l_n1(q);

            // -eᵀ [n1×] d²n0/(dxp dxq): d²n0_k/(dxp dxq) is in d2n0_local
            if (lp0 >= 0 && lq0 >= 0) {
                for (int k = 0; k < 3; k++) {
                    H_sin(p, q) +=
                        neg_e_cross_n1(k) * d2n0_local(3 * lp0 + k, lq0);
                }
            }
            // +eᵀ [n0×] d²n1/(dxp dxq)
            if (lp1 >= 0 && lq1 >= 0) {
                for (int k = 0; k < 3; k++) {
                    H_sin(p, q) += e_cross_n0(k) * d2n1_local(3 * lp1 + k, lq1);
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Assemble the Hessian of θ = atan2(s, c)
    //
    // Since s² + c² = 1 (n0, n1 are unit normals and e ⊥ plane spanned by
    // n0 and n1, so (n0×n1)·e = ||n0×n1|| = sin(φ):
    //
    // H_θ = c·H_s - s·H_c            (second-derivative terms)
    //      + ∇c·(∇s)ᵀ - ∇s·(∇c)ᵀ     (product-rule cross terms, antisymmetric)
    //      - 2·∇θ·(s·∇s + c·∇c)ᵀ     (denominator derivative)
    // -------------------------------------------------------------------------

    Matrix12 H_theta = cos_theta * H_sin - sin_theta * H_cos
        + dcos_dx * dsin_dx.transpose() - dsin_dx * dcos_dx.transpose()
        - dtheta_dx
            * ((T(2) * sin_theta) * dsin_dx + (T(2) * cos_theta) * dcos_dx)
                  .transpose();

    return H_theta;
}

// clang-format off
#define IPC_INSTANTIATE_ANGLE(T)                                              \
    template T dihedral_angle<T>(Eigen::ConstRef<Eigen::Vector3<T>>, Eigen::ConstRef<Eigen::Vector3<T>>, Eigen::ConstRef<Eigen::Vector3<T>>, Eigen::ConstRef<Eigen::Vector3<T>>); \
    template Eigen::Vector<T, 12> dihedral_angle_gradient<T>(Eigen::ConstRef<Eigen::Vector3<T>>, Eigen::ConstRef<Eigen::Vector3<T>>, Eigen::ConstRef<Eigen::Vector3<T>>, Eigen::ConstRef<Eigen::Vector3<T>>); \
    template Eigen::Matrix<T, 12, 12> dihedral_angle_hessian<T>(Eigen::ConstRef<Eigen::Vector3<T>>, Eigen::ConstRef<Eigen::Vector3<T>>, Eigen::ConstRef<Eigen::Vector3<T>>, Eigen::ConstRef<Eigen::Vector3<T>>)

IPC_INSTANTIATE_ANGLE(float);
IPC_INSTANTIATE_ANGLE(double);
#ifdef IPC_TOOLKIT_WITH_SIMD
IPC_INSTANTIATE_ANGLE(SimdBatch<float>);
IPC_INSTANTIATE_ANGLE(SimdBatch<double>);
#endif

#undef IPC_INSTANTIATE_ANGLE
// clang-format on

} // namespace ipc::detail
