#include "tangent_basis.hpp"

#include <ipc/utils/simd.hpp>

#include <Eigen/Geometry>

namespace ipc::detail {

// ========================================================================
// Point - Point

template <typename T, int dim>
Eigen::Matrix<T, dim*(dim - 1), 2 * dim> point_point_tangent_basis_jacobian(
    Eigen::ConstRef<Eigen::Vector<T, dim>> p0,
    Eigen::ConstRef<Eigen::Vector<T, dim>> p1)
{
    static_assert(
        dim == 2 || dim == 3, "point-point tangent basis is only 2D or 3D");

    Eigen::Matrix<T, dim*(dim - 1), 2 * dim> J;
    if constexpr (dim == 2) {
        autogen::point_point_tangent_basis_2D_jacobian<T>(
            p0[0], p0[1], p1[0], p1[1], J.data());
    } else {
        autogen::point_point_tangent_basis_3D_jacobian<T>(
            p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], J.data());
    }
    return J;
}

// ========================================================================
// Point - Edge

template <typename T, int dim>
Eigen::Matrix<T, dim*(dim - 1), 3 * dim> point_edge_tangent_basis_jacobian(
    Eigen::ConstRef<Eigen::Vector<T, dim>> p,
    Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
    Eigen::ConstRef<Eigen::Vector<T, dim>> e1)
{
    static_assert(
        dim == 2 || dim == 3, "point-edge tangent basis is only 2D or 3D");

    Eigen::Matrix<T, dim*(dim - 1), 3 * dim> J;
    if constexpr (dim == 2) {
        autogen::point_edge_tangent_basis_2D_jacobian<T>(
            p[0], p[1], e0[0], e0[1], e1[0], e1[1], J.data());
    } else {
        autogen::point_edge_tangent_basis_3D_jacobian<T>(
            p[0], p[1], p[2], e0[0], e0[1], e0[2], e1[0], e1[1], e1[2],
            J.data());
    }
    return J;
}

// ========================================================================
// Edge - Edge

template <typename T>
Eigen::Matrix<T, 6, 12> edge_edge_tangent_basis_jacobian(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    Eigen::Matrix<T, 6, 12> J;
    autogen::edge_edge_tangent_basis_jacobian<T>(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], J.data());
    return J;
}

// ========================================================================
// Point - Triangle

template <typename T>
Eigen::Matrix<T, 6, 12> point_triangle_tangent_basis_jacobian(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2)
{
    Eigen::Matrix<T, 6, 12> J;
    autogen::point_triangle_tangent_basis_jacobian<T>(
        p[0], p[1], p[2], t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0],
        t2[1], t2[2], J.data());
    return J;
}

// ========================================================================
// Explicit instantiations

#define IPC_INSTANTIATE_TANGENT_BASIS_ND(T, dim)                               \
    template Eigen::Matrix<T, dim, dim - 1> point_point_tangent_basis<T, dim>( \
        Eigen::ConstRef<Eigen::Vector<T, dim>>,                                \
        Eigen::ConstRef<Eigen::Vector<T, dim>>);                               \
    template Eigen::Matrix<T, dim*(dim - 1), 2 * dim>                          \
        point_point_tangent_basis_jacobian<T, dim>(                            \
            Eigen::ConstRef<Eigen::Vector<T, dim>>,                            \
            Eigen::ConstRef<Eigen::Vector<T, dim>>);                           \
    template Eigen::Matrix<T, dim, dim - 1> point_edge_tangent_basis<T, dim>(  \
        Eigen::ConstRef<Eigen::Vector<T, dim>>,                                \
        Eigen::ConstRef<Eigen::Vector<T, dim>>,                                \
        Eigen::ConstRef<Eigen::Vector<T, dim>>);                               \
    template Eigen::Matrix<T, dim*(dim - 1), 3 * dim>                          \
    point_edge_tangent_basis_jacobian<T, dim>(                                 \
        Eigen::ConstRef<Eigen::Vector<T, dim>>,                                \
        Eigen::ConstRef<Eigen::Vector<T, dim>>,                                \
        Eigen::ConstRef<Eigen::Vector<T, dim>>)

IPC_INSTANTIATE_TANGENT_BASIS_ND(float, 2);
IPC_INSTANTIATE_TANGENT_BASIS_ND(float, 3);
IPC_INSTANTIATE_TANGENT_BASIS_ND(double, 2);
IPC_INSTANTIATE_TANGENT_BASIS_ND(double, 3);
#ifdef IPC_TOOLKIT_WITH_SIMD
IPC_INSTANTIATE_TANGENT_BASIS_ND(SimdBatch<float>, 2);
IPC_INSTANTIATE_TANGENT_BASIS_ND(SimdBatch<float>, 3);
IPC_INSTANTIATE_TANGENT_BASIS_ND(SimdBatch<double>, 2);
IPC_INSTANTIATE_TANGENT_BASIS_ND(SimdBatch<double>, 3);
#endif

#undef IPC_INSTANTIATE_TANGENT_BASIS_ND

#define IPC_INSTANTIATE_TANGENT_BASIS_3D(T)                                    \
    template Eigen::Matrix<T, 3, 2> edge_edge_tangent_basis<T>(                \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>);                                   \
    template Eigen::Matrix<T, 6, 12> edge_edge_tangent_basis_jacobian<T>(      \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>);                                   \
    template Eigen::Matrix<T, 3, 2> point_triangle_tangent_basis<T>(           \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>);                                   \
    template Eigen::Matrix<T, 6, 12> point_triangle_tangent_basis_jacobian<T>( \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>,                                    \
        Eigen::ConstRef<Eigen::Vector3<T>>)

IPC_INSTANTIATE_TANGENT_BASIS_3D(float);
IPC_INSTANTIATE_TANGENT_BASIS_3D(double);
#ifdef IPC_TOOLKIT_WITH_SIMD
IPC_INSTANTIATE_TANGENT_BASIS_3D(SimdBatch<float>);
IPC_INSTANTIATE_TANGENT_BASIS_3D(SimdBatch<double>);
#endif

#undef IPC_INSTANTIATE_TANGENT_BASIS_3D

} // namespace ipc::detail

namespace ipc::autogen {

namespace {
    /// @brief Compute the power of 1.5 of a number.
    /// @param x Number to compute the power of 1.5
    /// @return x^(1.5)
    template <typename T> inline T pow_1_5(T x) { return x * ipc::sqrt(x); }
} // namespace

// J is (2×4) flattened in column-major order
template <typename T>
void point_point_tangent_basis_2D_jacobian(
    T p0_x, T p0_y, T p1_x, T p1_y, T J[8])
{
    const T t0 = p0_x - p1_x;
    const T t1 = p0_y - p1_y;
    const T t2 = t0 * t0;
    const T t3 = t1 * t1;
    const T t4 = t2 + t3;
    const T t5 = t0 * t1 / pow_1_5(t4);
    const T t6 = -t5;
    const T t7 = T(1) / ipc::sqrt(t4);
    const T t8 = T(1) / t4;
    const T t9 = t2 * t8 - 1;
    const T t10 = t3 * t8 - 1;
    J[0] = t6;
    J[1] = t7 * t9;
    J[2] = -t10 * t7;
    J[3] = t5;
    J[4] = t5;
    J[5] = -t7 * t9;
    J[6] = t10 * t7;
    J[7] = t6;
}

// J is (6×6) flattened in column-major order
template <typename T>
void point_point_tangent_basis_3D_jacobian(
    T p0_x, T p0_y, T p0_z, T p1_x, T p1_y, T p1_z, T J[36])
{
    const T t0 = p0_x - p1_x;
    const T t1 = t0 * t0;
    const T t2 = p0_z - p1_z;
    const T t3 = t2 * t2;
    const T t4 = p0_y - p1_y;
    const T t5 = t4 * t4;
    const auto t6 = (t1 + t3) < (t3 + t5);
    const T t7 = -t2;
    const T t8 = select(t6, T(0), t7);
    const T t9 =
        select(t6, T(0), t3) + select(t6, t3, T(0)) + select(t6, t5, t1);
    const T t10 = T(1) / pow_1_5(t9);
    const T t11 = T(0.5) * select(t6, T(0), T(2) * t0);
    const T t12 = t10 * t11;
    const T t13 = select(t6, t2, T(0));
    const T t14 = T(1) / ipc::sqrt(t9);
    const T t15 = select(t6, T(0), T(1));
    const T t16 = -t4;
    const T t17 = select(t6, t16, t0);
    const T t18 = T(1) / t9;
    const T t19 = t17 * t18;
    const T t20 = t0 * t13 - t4 * t8;
    const T t21 = -t20;
    const auto t22 = t1 + (t7 * t7) < (t16 * t16) + t3;
    const T t23 = select(t22, T(0), t7);
    const T t24 = -t0;
    const T t25 = select(t22, t16, t0);
    const T t26 = -t23 * t7 + t24 * t25;
    const T t27 = -t13 * t2 + t17 * t4;
    const T t28 = -t27;
    const T t29 = t21 * t21 + t26 * t26 + t28 * t28;
    const T t30 = T(1) / ipc::sqrt(t29);
    const T t31 = t15 * t4;
    const T t32 = T(1) / t29;
    const T t33 = select(t22, t2, T(0));
    const T t34 = -t16 * t23 + t24 * t33;
    const T t35 = t33 * t34;
    const T t36 = t16 * t25 - t33 * t7;
    const T t37 = select(t22, T(0), T(1));
    const T t38 = t16 * t37;
    const T t39 = t24 * t37;
    const T t40 = -t25;
    const T t41 = -t26;
    const T t42 = t32 * (-t35 + t36 * t38 + t41 * (-t39 - t40));
    const T t43 = T(0.5) * select(t6, T(2) * t4, T(0));
    const T t44 = t10 * t43;
    const T t45 = select(t6, T(-1), T(0));
    const T t46 = -t0 * t17 + t2 * t8;
    const T t47 = t20 * t20 + t27 * t27 + t46 * t46;
    const T t48 = T(1) / ipc::sqrt(t47);
    const T t49 = T(1) / t47;
    const T t50 = t0 * t45;
    const T t51 = t17 + t4 * t45;
    const T t52 = select(t22, T(-1), T(0));
    const T t53 = t24 * t52;
    const T t54 = t32 * (t23 * t34 + t36 * (t16 * t52 + t40) - t41 * t53);
    const T t55 = -t8;
    const T t56 = select(t6, T(0), T(-1));
    const T t57 = T(0.5) * t8;
    const T t58 = 2 * t2;
    const T t59 = t18 * (select(t22, T(0), t58) + select(t22, t58, T(0)));
    const T t60 = select(t6, T(1), T(0));
    const T t61 = T(0.5) * t13;
    const T t62 = T(0.5) * t10 * t17;
    const T t63 = select(t22, T(1), T(0));
    const T t64 = select(t22, T(0), T(-1));
    const T t65 = t16 * t64;
    const T t66 = -t33 + t63 * t7;
    const T t67 = t28 * t32;
    const T t68 = t0 * t60;
    const T t69 = t2 * t56 + t8;
    const T t70 = t21 * (t4 * t56 - t68) + t26 * t69 - t28 * t66;
    const T t71 = t4 * t56;
    const T t72 = select(t6, T(0), T(2) * t24);
    const T t73 = t10 * t72;
    const T t74 = T(0.5) * t19;
    const T t75 = t24 * t64 + t25;
    const T t76 = t35 + t36 * t65 - t41 * t75;
    const T t77 = select(t6, T(2) * t16, T(0));
    const T t78 = t10 * t77;
    const T t79 = -t17 + t4 * t60;
    const T t80 = t0 * t46 * t60 - t20 * t8 - t27 * t79;
    const T t81 = t20 * t49;
    const T t82 = 2 * t7;
    const T t83 = t18 * (select(t22, T(0), t82) + select(t22, t82, T(0)));
    const T t84 = t33 + t52 * t7;
    const T t85 = t15 * t2;
    const T t86 = t21 * (t15 * t4 - t50) - t26 * (-t55 - t85) - t28 * t84;
    J[0] = -t12 * t8;
    J[1] = -t12 * t13;
    J[2] = t14 * (-t11 * t19 + t15);
    J[3] = -t30 * (t28 * t42 + t31);
    J[4] = t30 * (t25 + t26 * t42 - t39);
    J[5] = -t30 * (t13 + t21 * t42);
    J[6] = -t44 * t8;
    J[7] = -t13 * t44;
    J[8] = t14 * (-t19 * t43 + t45);
    J[9] = t48 * (t27 * t49 * (t21 * t8 - t26 * t50 - t28 * t51) - t51);
    J[10] = t30 * (t26 * t54 + t50);
    J[11] = t30 * (-t21 * t54 - t55);
    J[12] = t14 * (t56 - t57 * t59);
    J[13] = t14 * (-t59 * t61 + t60);
    J[14] = -t62 * (select(t6, T(0), t58) + select(t6, t58, T(0)));
    J[15] = -t30
        * (t66
           + t67
               * (t34 * (t24 * t63 - t65) - t36 * t66
                  + t41 * (-t23 + t64 * t7)));
    J[16] = t48 * (t46 * t49 * t70 - t69);
    J[17] = t48 * (t20 * t49 * t70 - t68 + t71);
    J[18] = -t57 * t73;
    J[19] = -t61 * t73;
    J[20] = t14 * (t56 - t72 * t74);
    J[21] = -t30 * (t67 * t76 + t71);
    J[22] = t30 * (t26 * t32 * t76 - t75);
    J[23] = t30 * (t13 - t21 * t32 * t76);
    J[24] = -t57 * t78;
    J[25] = -t61 * t78;
    J[26] = t14 * (t60 - t74 * t77);
    J[27] = -t48 * (t27 * t49 * t80 + t79);
    J[28] = t30 * (-t26 * t32 * (t21 * t8 + t26 * t68 + t28 * t79) + t68);
    J[29] = -t48 * (t8 + t80 * t81);
    J[30] = t14 * (t15 - t57 * t83);
    J[31] = t14 * (t45 - t61 * t83);
    J[32] = -t62 * (select(t6, T(0), t82) + select(t6, t82, T(0)));
    J[33] = -t30
        * (t67 * (t34 * (-t38 + t53) - t36 * t84 + t41 * (t23 + t37 * t7))
           + t84);
    J[34] = t48 * (t46 * t49 * t86 + t8 - t85);
    J[35] = t48 * (t31 - t50 + t81 * t86);
}

// J is (2×6) flattened in column-major order
template <typename T>
void point_edge_tangent_basis_2D_jacobian(
    T p_x, T p_y, T e0_x, T e0_y, T e1_x, T e1_y, T J[12])
{
    const T t0 = e0_x - e1_x;
    const T t1 = t0 * t0;
    const T t2 = e0_y - e1_y;
    const T t3 = t2 * t2;
    const T t4 = t1 + t3;
    const T t5 = T(1) / ipc::sqrt(t4);
    const T t6 = T(1) / t4;
    const T t7 = t1 * t6 - 1;
    const T t8 = t0 * t2 / (t4 * ipc::sqrt(t4));
    const T t9 = t3 * t6 - 1;
    const T t10 = -t8;
    J[0] = 0;
    J[1] = 0;
    J[2] = 0;
    J[3] = 0;
    J[4] = t5 * t7;
    J[5] = t8;
    J[6] = t8;
    J[7] = t5 * t9;
    J[8] = -t5 * t7;
    J[9] = t10;
    J[10] = t10;
    J[11] = -t5 * t9;
}

// J is (6×9) flattened in column-major order
template <typename T>
void point_edge_tangent_basis_3D_jacobian(
    T p_x,
    T p_y,
    T p_z,
    T e0_x,
    T e0_y,
    T e0_z,
    T e1_x,
    T e1_y,
    T e1_z,
    T J[54])
{
    const T t0 = -e1_y;
    const T t1 = e0_y + t0;
    const T t2 = -e1_x;
    const T t3 = e0_x + t2;
    const T t4 = -p_y;
    const T t5 = e0_y + t4;
    const T t6 = -p_x;
    const T t7 = e0_x + t6;
    const T t8 = -t1 * t7 + t3 * t5;
    const T t9 = -e1_z;
    const T t10 = e0_z + t9;
    const T t11 = -t3;
    const T t12 = -p_z;
    const T t13 = e0_z + t12;
    const T t14 = -t13;
    const T t15 = -t7;
    const T t16 = -t10;
    const T t17 = t11 * t14 - t15 * t16;
    const T t18 = t1 * t8 + t10 * t17;
    const T t19 = t1 * t13 - t10 * t5;
    const T t20 = -t10 * t7 + t13 * t3;
    const T t21 = t19 * t19 + t8 * t8;
    const T t22 = t20 * t20 + t21;
    const T t23 = T(1) / pow_1_5(t22);
    const T t24 = t19 * t23;
    const T t25 = t17 * t17 + t21;
    const T t26 = T(1) / ipc::sqrt(t25);
    const T t27 = -t5;
    const T t28 = -t1;
    const T t29 = t11 * t27 - t15 * t28;
    const T t30 = -t17;
    const T t31 = T(1) / t25;
    const T t32 = t17 * t31;
    const T t33 = -e0_z;
    const T t34 = e1_z + t33;
    const T t35 = T(1) / ipc::sqrt(t22);
    const T t36 = -e0_y;
    const T t37 = T(1) / t22;
    const T t38 = t37 * t8;
    const T t39 = t10 * t19 - t3 * t8;
    const T t40 = t19 * t37;
    const T t41 = t20 * t23;
    const T t42 = t14 * t28 - t16 * t27;
    const T t43 = t31 * t8;
    const T t44 = t19 * t31;
    const T t45 = -e0_x;
    const T t46 = t1 * t19 + t17 * t3;
    const T t47 = t20 * t37;
    const T t48 = t23 * t8;
    const T t49 = t3 * t3;
    const T t50 = t1 * t1;
    const T t51 = t10 * t10;
    const T t52 = t49 + t50 + t51;
    const T t53 = T(1) / ipc::sqrt(t52);
    const T t54 = T(1) / t52;
    const T t55 = t49 * t54 - 1;
    const T t56 = T(1) / pow_1_5(t52);
    const T t57 = t3 * t56;
    const T t58 = t1 * t57;
    const T t59 = t10 * t57;
    const T t60 = e1_y + t4;
    const T t61 = e1_z + t12;
    const T t62 = t17 * t61 + t60 * t8;
    const T t63 = p_z + t9;
    const T t64 = t50 * t54 - 1;
    const T t65 = t1 * t10 * t56;
    const T t66 = e1_x + t6;
    const T t67 = t19 * t61 - t66 * t8;
    const T t68 = t51 * t54 - 1;
    const T t69 = t17 * t66 + t19 * t60;
    const T t70 = -t58;
    const T t71 = -t59;
    const T t72 = t13 * t17 + t5 * t8;
    const T t73 = -t65;
    const T t74 = -t13 * t19 + t7 * t8;
    const T t75 = p_x + t45;
    const T t76 = t17 * t7 + t19 * t5;
    J[0] = 0;
    J[1] = 0;
    J[2] = 0;
    J[3] = -t18 * t24;
    J[4] = t26 * (t32 * (t1 * t29 + t16 * t30) + t34);
    J[5] = t35 * (-e1_y - t18 * t38 - t36);
    J[6] = 0;
    J[7] = 0;
    J[8] = 0;
    J[9] = t35 * (-t34 - t39 * t40);
    J[10] = t39 * t41;
    J[11] = -t26 * (t3 + t43 * (t10 * t42 + t11 * t29));
    J[12] = 0;
    J[13] = 0;
    J[14] = 0;
    J[15] = -t26 * (t1 + t44 * (t28 * t42 + t3 * t30));
    J[16] = t35 * (-e1_x - t45 - t46 * t47);
    J[17] = t46 * t48;
    J[18] = t53 * t55;
    J[19] = t58;
    J[20] = t59;
    J[21] = -t24 * t62;
    J[22] = t26 * (t32 * (t29 * t60 - t30 * t61) + t63);
    J[23] = t35 * (-p_y - t0 - t38 * t62);
    J[24] = t58;
    J[25] = t53 * t64;
    J[26] = t65;
    J[27] = t35 * (-t40 * t67 - t63);
    J[28] = t41 * t67;
    J[29] = -t26 * (t43 * (-t29 * t66 + t42 * t61) + t66);
    J[30] = t59;
    J[31] = t65;
    J[32] = t53 * t68;
    J[33] = -t26 * (t44 * (t30 * t66 - t42 * t60) + t60);
    J[34] = t35 * (-p_x - t2 - t47 * t69);
    J[35] = t48 * t69;
    J[36] = -t53 * t55;
    J[37] = t70;
    J[38] = t71;
    J[39] = t24 * t72;
    J[40] = t35 * (-p_z - t33 - t47 * t72);
    J[41] = -t26 * (t43 * (t13 * t30 + t27 * t29) + t5);
    J[42] = t70;
    J[43] = -t53 * t64;
    J[44] = t73;
    J[45] = -t26 * (t13 + t44 * (t14 * t42 + t29 * t7));
    J[46] = t41 * t74;
    J[47] = t35 * (-t38 * t74 - t75);
    J[48] = t71;
    J[49] = t73;
    J[50] = -t53 * t68;
    J[51] = t35 * (-p_y - t36 - t40 * t76);
    J[52] = t26 * (t32 * (t15 * t30 + t42 * t5) + t75);
    J[53] = -t48 * t76;
}

// J is (6×12) flattened in column-major order
template <typename T>
void edge_edge_tangent_basis_jacobian(
    T ea0_x,
    T ea0_y,
    T ea0_z,
    T ea1_x,
    T ea1_y,
    T ea1_z,
    T eb0_x,
    T eb0_y,
    T eb0_z,
    T eb1_x,
    T eb1_y,
    T eb1_z,
    T J[72])
{
    const T t0 = ea0_y - ea1_y;
    const T t1 = t0 * t0;
    const T t2 = ea0_x - ea1_x;
    const T t3 = t2 * t2;
    const T t4 = ea0_z - ea1_z;
    const T t5 = t4 * t4;
    const T t6 = t3 + t5;
    const T t7 = t1 + t6;
    const T t8 = T(1) / ipc::sqrt(t7);
    const T t9 = T(1) / t7;
    const T t10 = t3 * t9 - 1;
    const T t11 = T(1) / pow_1_5(t7);
    const T t12 = t0 * t2;
    const T t13 = t11 * t12;
    const T t14 = t2 * t4;
    const T t15 = t11 * t14;
    const T t16 = -t2;
    const T t17 = eb0_z - eb1_z;
    const T t18 = -t17;
    const T t19 = t16 * t18;
    const T t20 = -t4;
    const T t21 = eb0_x - eb1_x;
    const T t22 = -t21;
    const T t23 = t20 * t22;
    const T t24 = -t23;
    const T t25 = t19 + t24;
    const T t26 = -t25;
    const T t27 = -t0;
    const T t28 = t18 * t27;
    const T t29 = eb0_y - eb1_y;
    const T t30 = -t29;
    const T t31 = t20 * t30;
    const T t32 = t28 - t31;
    const T t33 = t16 * t26 - t27 * t32;
    const T t34 = t2 * t29;
    const T t35 = t0 * t21;
    const T t36 = -t35;
    const T t37 = t34 + t36;
    const T t38 = t0 * t17;
    const T t39 = t29 * t4;
    const T t40 = -t39;
    const T t41 = t38 + t40;
    const T t42 = t2 * t37 - t4 * t41;
    const T t43 = -t42;
    const T t44 = t16 * t30;
    const T t45 = t22 * t27;
    const T t46 = -t45;
    const T t47 = t44 + t46;
    const T t48 = -t20 * t26 + t27 * t47;
    const T t49 = t33 * t33 + t43 * t43 + t48 * t48;
    const T t50 = T(1) / ipc::sqrt(t49);
    const T t51 = t27 * t29;
    const T t52 = -t51;
    const T t53 = T(1) / t49;
    const T t54 = 2 * t19;
    const T t55 = t24 + t54;
    const T t56 = -t33;
    const T t57 = -t18 * t20 + t51;
    const T t58 = -t48;
    const T t59 = t16 * t47 - t20 * t32;
    const T t60 =
        t53 * (-t55 * t56 - t57 * t58 + t59 * (t16 * t29 - t44 + t45));
    const T t61 = t0 * t41 + t2 * t25;
    const T t62 = t0 * t37 + t25 * t4;
    const T t63 = t42 * t42 + t61 * t61 + t62 * t62;
    const T t64 = T(1) / ipc::sqrt(t63);
    const T t65 = 2 * t34;
    const T t66 = t36 + t65;
    const T t67 = T(1) / t63;
    const T t68 = t42 * t67;
    const T t69 = t1 * t9 - 1;
    const T t70 = t0 * t4;
    const T t71 = t11 * t70;
    const T t72 = 2 * t35;
    const T t73 = 2 * t38;
    const T t74 = t40 + t73;
    const T t75 = t17 * t4 + t2 * t21;
    const T t76 = t43 * t75;
    const T t77 = t34 - t72;
    const T t78 = t33 * t74 - t48 * t77 + t76;
    const T t79 = t67 * t78;
    const T t80 = t5 * t9 - 1;
    const T t81 = t21 * t4;
    const T t82 = 2 * t81;
    const T t83 = t16 * t21;
    const T t84 = t27 * t30;
    const T t85 = -t84;
    const T t86 = t83 + t85;
    const T t87 = 2 * t39;
    const T t88 = t17 * t2;
    const T t89 = t62 * t67;
    const T t90 =
        t53 * (-t56 * t86 + t58 * (t20 * t21 + t25) + t59 * (t28 - 2 * t31));
    const T t91 = -t13;
    const T t92 = -t15;
    const T t93 = -t17 * t20 + t84;
    const T t94 = -t19;
    const T t95 = t16 * t17 + t23 + t94;
    const T t96 = t48 * t53;
    const T t97 = t33 * t95 + t43 * t66 + t48 * t93;
    const T t98 = t61 * t67;
    const T t99 = -t71;
    const T t100 = -t33 * t74 + t48 * (t21 * t27 + t47) - t76;
    const T t101 = t16 * t22;
    const T t102 = t20 * t29 + t32;
    const T t103 =
        -t102 * t59 + t56 * (-t101 - t52) + t58 * (-t19 + 2 * t20 * t22);
    const T t104 = t103 * t53;
    const T t105 = t27 * t27;
    const T t106 = -t20 * t4;
    const T t107 = t105 + t106;
    const T t108 = t107 * t48 + t12 * t43 - t14 * t33;
    const T t109 = t16 * t59;
    const T t110 = t16 * t56;
    const T t111 = t33 * t53;
    const T t112 = t27 * t56;
    const T t113 = t27 * t58;
    const T t114 = t20 * t20;
    const T t115 = -t114;
    const T t116 = t16 * t2;
    const T t117 = t112 * t20 - t113 * t2 + t59 * (t115 + t116);
    const T t118 = t43 * t53;
    const T t119 = t16 * t16;
    const T t120 = t0 * t27;
    const T t121 = -t120;
    const T t122 = t119 + t121;
    const T t123 = t122 * t33 - t14 * t48 + t43 * t70;
    const T t124 = t20 * t58;
    const T t125 = t20 * t59;
    const T t126 = t0 * t109 - t110 * t20 + t58 * (-t115 - t120);
    const T t127 = t112 * t4 - t113 * t16 + t59 * (t106 + t119);
    const T t128 = t124 * t2 - t125 * t27 + t56 * (t105 - t116);
    J[0] = t10 * t8;
    J[1] = t13;
    J[2] = t15;
    J[3] = t50 * (t18 * t20 + t48 * t60 + t52);
    J[4] = t64 * (t35 - t65 + t68 * (t33 * t55 - t43 * t66 + t48 * t57));
    J[5] = t50 * (t23 + t33 * t60 - t54);
    J[6] = t13;
    J[7] = t69 * t8;
    J[8] = t71;
    J[9] = t64 * (t2 * t29 - t62 * t79 - t72);
    J[10] = t64 * (t68 * t78 + t75);
    J[11] = t64 * (t39 + t61 * t79 - t73);
    J[12] = t15;
    J[13] = t71;
    J[14] = t8 * t80;
    J[15] = t64
        * (t17 * t2 - t82
           - t89 * (t33 * t86 + t43 * (t38 - t87) - t48 * (-t82 + t88)));
    J[16] = t50 * (t0 * t17 - t43 * t90 - t87);
    J[17] = t50 * (t33 * t90 - t83 + t84);
    J[18] = -t10 * t8;
    J[19] = t91;
    J[20] = t92;
    J[21] = t50
        * (t17 * t20 + t85
           + t96 * (-t56 * t95 - t58 * t93 + t59 * (2 * t44 + t46)));
    J[22] = t64 * (t66 + t68 * t97);
    J[23] = t64 * (-t81 + 2 * t88 + t97 * t98);
    J[24] = t91;
    J[25] = -t69 * t8;
    J[26] = t99;
    J[27] = -t64 * (t100 * t89 + t77);
    J[28] = t64 * (t100 * t42 * t67 - t75);
    J[29] = t64 * (t100 * t98 + t74);
    J[30] = t92;
    J[31] = t99;
    J[32] = -t8 * t80;
    J[33] = t50 * (t103 * t96 + 2 * t23 + t94);
    J[34] = -t50 * (t102 + t104 * t43);
    J[35] = t50 * (-t101 + t104 * t33 + t51);
    J[36] = 0;
    J[37] = 0;
    J[38] = 0;
    J[39] = -t64 * (t1 + t108 * t89 + t5);
    J[40] = t64 * (t108 * t68 + t12);
    J[41] = t50 * (t111 * (-t107 * t58 + t109 * t27 - t110 * t4) + t14);
    J[42] = 0;
    J[43] = 0;
    J[44] = 0;
    J[45] = t50 * (t117 * t96 + t12);
    J[46] = -t50 * (t117 * t118 + t6);
    J[47] = t64 * (t70 - t98 * (t12 * t48 + t33 * t70 + t43 * t6));
    J[48] = 0;
    J[49] = 0;
    J[50] = 0;
    J[51] = t64 * (-t123 * t89 + t14);
    J[52] = t50 * (-t118 * (-t0 * t125 - t122 * t56 + t124 * t16) + t70);
    J[53] = t64 * (-t1 + t123 * t61 * t67 - t3);
    J[54] = 0;
    J[55] = 0;
    J[56] = 0;
    J[57] = t50 * (t114 + t121 + t126 * t96);
    J[58] = -t50 * (t118 * t126 + t12);
    J[59] = t50 * (t126 * t33 * t53 - t14);
    J[60] = 0;
    J[61] = 0;
    J[62] = 0;
    J[63] = t50 * (-t12 + t127 * t48 * t53);
    J[64] = t50 * (-t118 * t127 + t6);
    J[65] = t50 * (t127 * t33 * t53 - t70);
    J[66] = 0;
    J[67] = 0;
    J[68] = 0;
    J[69] = t50 * (t128 * t48 * t53 - t14);
    J[70] = -t50 * (t118 * t128 + t70);
    J[71] = t50 * (t105 + t111 * t128 - t116);
}

// J is (6×12) flattened in column-major order
template <typename T>
void point_triangle_tangent_basis_jacobian(
    T p_x,
    T p_y,
    T p_z,
    T t0_x,
    T t0_y,
    T t0_z,
    T t1_x,
    T t1_y,
    T t1_z,
    T t2_x,
    T t2_y,
    T t2_z,
    T J[72])
{
    const T t0 = t0_y - t1_y;
    const T t1 = t0 * t0;
    const T t2 = t0_x - t1_x;
    const T t3 = t2 * t2;
    const T t4 = t0_z - t1_z;
    const T t5 = t4 * t4;
    const T t6 = t3 + t5;
    const T t7 = t1 + t6;
    const T t8 = (T(1) / ipc::sqrt(t7));
    const T t9 = T(1) / t7;
    const T t10 = t3 * t9 - 1;
    const T t11 = T(1) / pow_1_5(t7);
    const T t12 = t0 * t2;
    const T t13 = t11 * t12;
    const T t14 = t2 * t4;
    const T t15 = t11 * t14;
    const T t16 = -t2;
    const T t17 = -t2_z;
    const T t18 = t0_z + t17;
    const T t19 = -t18;
    const T t20 = t16 * t19;
    const T t21 = -t2_x;
    const T t22 = t0_x + t21;
    const T t23 = -t22;
    const T t24 = -t4;
    const T t25 = t23 * t24;
    const T t26 = t20 - t25;
    const T t27 = -t26;
    const T t28 = -t0;
    const T t29 = t19 * t28;
    const T t30 = -t2_y;
    const T t31 = t0_y + t30;
    const T t32 = -t31;
    const T t33 = t24 * t32;
    const T t34 = t29 - t33;
    const T t35 = t16 * t27 - t28 * t34;
    const T t36 = t2 * t31;
    const T t37 = t0 * t22;
    const T t38 = -t37;
    const T t39 = t36 + t38;
    const T t40 = t0 * t18;
    const T t41 = -t31 * t4;
    const T t42 = t40 + t41;
    const T t43 = t2 * t39 - t4 * t42;
    const T t44 = -t43;
    const T t45 = t16 * t32;
    const T t46 = t23 * t28;
    const T t47 = -t46;
    const T t48 = t45 + t47;
    const T t49 = -t24 * t27 + t28 * t48;
    const T t50 = t35 * t35 + t44 * t44 + t49 * t49;
    const T t51 = T(1) / ipc::sqrt(t50);
    const T t52 = t17 + t1_z;
    const T t53 = -t52;
    const T t54 = t1_y + t30;
    const T t55 = t28 * t54;
    const T t56 = T(1) / t50;
    const T t57 = -t24 * t53 + t55;
    const T t58 = -t49;
    const T t59 = -t45 + t46;
    const T t60 = t16 * t48 - t24 * t34;
    const T t61 = t16 * t53 + t26;
    const T t62 = -t35;
    const T t63 = t56 * (-t57 * t58 + t60 * (t16 * t54 + t59) - t61 * t62);
    const T t64 = t2 * t54 + t39;
    const T t65 = t0 * t42 + t2 * t26;
    const T t66 = t0 * t39 + t26 * t4;
    const T t67 = t43 * t43 + t65 * t65 + t66 * t66;
    const T t68 = T(1) / ipc::sqrt(t67);
    const T t69 = t18 * t2;
    const T t70 = t22 * t4;
    const T t71 = -t70;
    const T t72 = T(1) / t67;
    const T t73 = t1 * t9 - 1;
    const T t74 = t0 * t4;
    const T t75 = t11 * t74;
    const T t76 = t1_x + t21;
    const T t77 = t2 * t76 + t4 * t52;
    const T t78 = t0 * t52 + t42;
    const T t79 = t35 * t78 + t44 * t77 + t49 * (-t28 * t76 + t59);
    const T t80 = t72 * t79;
    const T t81 = t5 * t9 - 1;
    const T t82 = t16 * t76;
    const T t83 = -t54;
    const T t84 = t28 * t83;
    const T t85 = t82 - t84;
    const T t86 = t24 * t83 - t29 + t33;
    const T t87 = t4 * t76 - t69 + t70;
    const T t88 = t66 * t72;
    const T t89 = t56 * (t58 * (t24 * t76 + t26) - t60 * t86 - t62 * t85);
    const T t90 = -t13;
    const T t91 = -t15;
    const T t92 = t28 * t32;
    const T t93 = -t18 * t24 + t92;
    const T t94 = -t20;
    const T t95 = t16 * t18 + t25 + t94;
    const T t96 = t49 * t56;
    const T t97 = 2 * t36 + t38;
    const T t98 = t35 * t95 + t44 * t97 + t49 * t93;
    const T t99 = t65 * t72;
    const T t100 = -t75;
    const T t101 = t18 * t4 + t2 * t22;
    const T t102 = 2 * t40 + t41;
    const T t103 = -t101 * t44 - t102 * t35 + t49 * (t22 * t28 + t48);
    const T t104 = t16 * t23;
    const T t105 = t24 * t31 + t34;
    const T t106 =
        -t105 * t60 + t58 * (-t20 + 2 * t23 * t24) + t62 * (-t104 + t28 * t31);
    const T t107 = t106 * t56;
    const T t108 = t24 * t24;
    const T t109 = t0 * t28;
    const T t110 = t0 * t16 * t60 - t16 * t24 * t62 + t58 * (t108 - t109);
    const T t111 =
        -t16 * t28 * t58 + t28 * t4 * t62 + t60 * (t16 * t16 - t24 * t4);
    const T t112 = t28 * t28;
    const T t113 = t16 * t2;
    const T t114 = t2 * t24 * t58 - t24 * t28 * t60 + t62 * (t112 - t113);
    const T t115 = t114 * t56;
    J[0] = 0;
    J[1] = 0;
    J[2] = 0;
    J[3] = 0;
    J[4] = 0;
    J[5] = 0;
    J[6] = 0;
    J[7] = 0;
    J[8] = 0;
    J[9] = 0;
    J[10] = 0;
    J[11] = 0;
    J[12] = 0;
    J[13] = 0;
    J[14] = 0;
    J[15] = 0;
    J[16] = 0;
    J[17] = 0;
    J[18] = t10 * t8;
    J[19] = t13;
    J[20] = t15;
    J[21] = t51 * (t24 * t53 + t49 * t63 - t55);
    J[22] = t51 * (-t44 * t63 - t64);
    J[23] = t68
        * (-t2 * t52 + t65 * t72 * (t35 * t61 - t44 * t64 + t49 * t57) - t69
           - t71);
    J[24] = t13;
    J[25] = t73 * t8;
    J[26] = t75;
    J[27] = -t68 * (t0 * t76 - t36 + t37 + t66 * t80);
    J[28] = t68 * (t43 * t80 + t77);
    J[29] = t68 * (t65 * t72 * t79 - t78);
    J[30] = t15;
    J[31] = t75;
    J[32] = t8 * t81;
    J[33] = -t68 * (t87 + t88 * (t35 * t85 - t44 * t86 + t49 * t87));
    J[34] = -t51 * (t44 * t89 + t86);
    J[35] = t51 * (t35 * t89 - t82 + t84);
    J[36] = -t10 * t8;
    J[37] = t90;
    J[38] = t91;
    J[39] = t51
        * (t18 * t24 - t92
           + t96 * (-t58 * t93 + t60 * (2 * t45 + t47) - t62 * t95));
    J[40] = t68 * (t43 * t72 * t98 + t97);
    J[41] = t68 * (2 * t69 + t71 + t98 * t99);
    J[42] = t90;
    J[43] = -t73 * t8;
    J[44] = t100;
    J[45] = -t68 * (t103 * t88 + t36 - 2 * t37);
    J[46] = t68 * (-t101 + t103 * t43 * t72);
    J[47] = t68 * (t102 + t103 * t99);
    J[48] = t91;
    J[49] = t100;
    J[50] = -t8 * t81;
    J[51] = t51 * (t106 * t96 + 2 * t25 + t94);
    J[52] = -t51 * (t105 + t107 * t44);
    J[53] = t51 * (-t104 + t107 * t35 + t28 * t31);
    J[54] = 0;
    J[55] = 0;
    J[56] = 0;
    J[57] = t51 * (t108 - t109 + t110 * t96);
    J[58] = -t51 * (t110 * t44 * t56 + t12);
    J[59] = t51 * (t110 * t35 * t56 - t14);
    J[60] = 0;
    J[61] = 0;
    J[62] = 0;
    J[63] = t51 * (t111 * t49 * t56 - t12);
    J[64] = t51 * (-t111 * t44 * t56 + t6);
    J[65] = t51 * (t111 * t35 * t56 - t74);
    J[66] = 0;
    J[67] = 0;
    J[68] = 0;
    J[69] = t51 * (t114 * t49 * t56 - t14);
    J[70] = -t51 * (t115 * t44 + t74);
    J[71] = t51 * (t112 - t113 + t115 * t35);
}

#define IPC_INSTANTIATE_TANGENT_BASIS_AUTOGEN(T)                               \
    template void point_point_tangent_basis_2D_jacobian<T>(T, T, T, T, T[8]);  \
    template void point_point_tangent_basis_3D_jacobian<T>(                    \
        T, T, T, T, T, T, T[36]);                                              \
    template void point_edge_tangent_basis_2D_jacobian<T>(                     \
        T, T, T, T, T, T, T[12]);                                              \
    template void point_edge_tangent_basis_3D_jacobian<T>(                     \
        T, T, T, T, T, T, T, T, T, T[54]);                                     \
    template void edge_edge_tangent_basis_jacobian<T>(                         \
        T, T, T, T, T, T, T, T, T, T, T, T, T[72]);                            \
    template void point_triangle_tangent_basis_jacobian<T>(                    \
        T, T, T, T, T, T, T, T, T, T, T, T, T[72])

IPC_INSTANTIATE_TANGENT_BASIS_AUTOGEN(float);
IPC_INSTANTIATE_TANGENT_BASIS_AUTOGEN(double);
#ifdef IPC_TOOLKIT_WITH_SIMD
IPC_INSTANTIATE_TANGENT_BASIS_AUTOGEN(SimdBatch<float>);
IPC_INSTANTIATE_TANGENT_BASIS_AUTOGEN(SimdBatch<double>);
#endif

#undef IPC_INSTANTIATE_TANGENT_BASIS_AUTOGEN

} // namespace ipc::autogen
