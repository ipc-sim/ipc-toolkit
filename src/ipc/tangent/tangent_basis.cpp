#include "tangent_basis.hpp"

#include <Eigen/Geometry>

namespace ipc {

namespace {
    /// @brief Compute the power of 1.5 of a number.
    /// @param x Number to compute the power of 1.5
    /// @return x^(1.5)
    inline double pow_1_5(double x) { return x * std::sqrt(x); }
} // namespace

// ============================================================================
// Point - Point

MatrixMax<double, 3, 2> point_point_tangent_basis(
    Eigen::ConstRef<VectorMax3d> p0, Eigen::ConstRef<VectorMax3d> p1)
{
    const int dim = p0.size();
    assert(dim == p1.size());

    if (dim == 2) {
        const Eigen::Vector2d p0_to_p1 = (p1 - p0).normalized();
        return Eigen::Vector2d(-p0_to_p1.y(), p0_to_p1.x());
    } else {
        assert(dim == 3);

        const Eigen::Vector3d p0_to_p1 = p1 - p0;

        const Eigen::Vector3d cross_x =
            Eigen::Vector3d::UnitX().cross(p0_to_p1);
        const Eigen::Vector3d cross_y =
            Eigen::Vector3d::UnitY().cross(p0_to_p1);

        Eigen::Matrix<double, 3, 2> basis;
        if (cross_x.squaredNorm() > cross_y.squaredNorm()) {
            basis.col(0) = cross_x.normalized();
            basis.col(1) = p0_to_p1.cross(cross_x).normalized();
        } else {
            basis.col(0) = cross_y.normalized();
            basis.col(1) = p0_to_p1.cross(cross_y).normalized();
        }

        return basis;
    }
}

MatrixMax<double, 18, 2> point_point_tangent_basis_jacobian(
    Eigen::ConstRef<VectorMax3d> p0, Eigen::ConstRef<VectorMax3d> p1)
{
    const int dim = p0.size();
    assert(dim == p1.size());

    MatrixMax<double, 18, 2> J;
    if (dim == 2) {
        J.resize(8, 1);
        autogen::point_point_tangent_basis_2D_jacobian(
            p0[0], p0[1], p1[0], p1[1], J.data());
    } else {
        assert(dim == 3);
        J.resize(18, 2);
        autogen::point_point_tangent_basis_3D_jacobian(
            p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], J.data());
    }
    return J;
}

// ============================================================================
// Point - Edge

MatrixMax<double, 3, 2> point_edge_tangent_basis(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1)
{
    const int dim = p.size();
    assert(dim == e0.size() && dim == e1.size());

    if (dim == 2) {
        return (e1 - e0).normalized();
    } else {
        assert(dim == 3);

        const Eigen::Vector3d e = e1 - e0;

        Eigen::Matrix<double, 3, 2> basis;
        basis.col(0) = e.normalized();
        basis.col(1) = e.cross(Eigen::Vector3d(p - e0)).normalized();
        return basis;
    }
}

MatrixMax<double, 27, 2> point_edge_tangent_basis_jacobian(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1)
{
    const int dim = p.size();
    assert(dim == e0.size() && dim == e1.size());

    MatrixMax<double, 27, 2> J;
    if (dim == 2) {
        J.resize(12, 1);
        autogen::point_edge_tangent_basis_2D_jacobian(
            p[0], p[1], e0[0], e0[1], e1[0], e1[1], J.data());
    } else {
        assert(dim == 3);
        J.resize(27, 2);
        autogen::point_edge_tangent_basis_3D_jacobian(
            p[0], p[1], p[2], e0[0], e0[1], e0[2], e1[0], e1[1], e1[2],
            J.data());
    }
    return J;
}

// ============================================================================
// Edge - Edge

Eigen::Matrix<double, 3, 2> edge_edge_tangent_basis(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    const Eigen::Vector3d ea = ea1 - ea0; // Edge A direction
    const Eigen::Vector3d normal = ea.cross(eb1 - eb0);
    // The normal will be zero if the edges are parallel (i.e. coplanar).
    assert(normal.norm() != 0);

    Eigen::Matrix<double, 3, 2> basis;
    // The first basis vector is along edge A.
    basis.col(0) = ea.normalized();
    // The second basis vector is orthogonal to the first and the edge-edge
    // normal.
    basis.col(1) = normal.cross(ea).normalized();
    return basis;
}

Eigen::Matrix<double, 36, 2> edge_edge_tangent_basis_jacobian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    Eigen::Matrix<double, 36, 2> J;
    autogen::edge_edge_tangent_basis_jacobian(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], J.data());
    return J;
}

// ============================================================================
// Point - Triangle

Eigen::Matrix<double, 3, 2> point_triangle_tangent_basis(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2)
{
    const Eigen::Vector3d e0 = t1 - t0;
    const Eigen::Vector3d normal = e0.cross(t2 - t0);
    assert(normal.norm() != 0);

    Eigen::Matrix<double, 3, 2> basis;

    // The first basis vector is along first edge of the triangle.
    basis.col(0) = e0.normalized();
    // The second basis vector is orthogonal to the first and the triangle
    // normal.
    basis.col(1) = normal.cross(e0).normalized();

    return basis;
}

Eigen::Matrix<double, 36, 2> point_triangle_tangent_basis_jacobian(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2)
{
    Eigen::Matrix<double, 36, 2> J;
    autogen::point_triangle_tangent_basis_jacobian(
        p[0], p[1], p[2], t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0],
        t2[1], t2[2], J.data());
    return J;
}

// ============================================================================

namespace autogen {

    // J is (8×1) flattened in column-major order
    void point_point_tangent_basis_2D_jacobian(
        double p0_x, double p0_y, double p1_x, double p1_y, double J[8])
    {
        const double t0 = p0_x - p1_x;
        const double t1 = p0_y - p1_y;
        const double t2 = t0 * t0;
        const double t3 = t1 * t1;
        const double t4 = t2 + t3;
        const double t5 = t0 * t1 / pow_1_5(t4);
        const double t6 = -t5;
        const double t7 = 1.0 / std::sqrt(t4);
        const double t8 = 1.0 / t4;
        const double t9 = t2 * t8;
        const double t10 = t3 * t8;
        J[0] = t6;
        J[1] = t7 * (t9 - 1);
        J[2] = t7 * (1 - t10);
        J[3] = t5;
        J[4] = t5;
        J[5] = t7 * (1 - t9);
        J[6] = t7 * (t10 - 1);
        J[7] = t6;
    }

    // J is (18×2) flattened in column-major order
    void point_point_tangent_basis_3D_jacobian(
        double p0_x,
        double p0_y,
        double p0_z,
        double p1_x,
        double p1_y,
        double p1_z,
        double J[36])
    {
        const double t0 = p0_x - p1_x;
        const double t1 = t0 * t0;
        const double t2 = p0_z - p1_z;
        const double t3 = t2 * t2;
        const double t4 = p0_y - p1_y;
        const double t5 = t4 * t4;
        const bool t6 = (t1 + t3) < (t3 + t5);
        const double t7 = -p0_z + p1_z;
        const double t8 = t6 ? 0 : t7;
        const double t9 = t3 + (t6 ? t5 : t1);
        const double t10 = 1.0 / pow_1_5(t9);
        const double t11 = t6 ? 0 : t0;
        const double t12 = t10 * t11;
        const double t13 = t6 ? t2 : 0;
        const double t14 = 1 / std::sqrt(t9);
        const double t15 = int(!t6);
        const double t16 = -p0_y + p1_y;
        const double t17 = t6 ? t16 : t0;
        const double t18 = 1.0 / t9;
        const double t19 = t17 * t18;
        const double t20 = t6 ? t4 : 0;
        const double t21 = t10 * t20;
        const double t22 = -int(t6);
        const double t23 = -int(!t6);
        const double t24 = 0.5 * t8;
        const bool t25 = t1 + (t7 * t7) < (t16 * t16) + t3;
        const double t26 = 2 * t2;
        const double t27 = t18 * t26;
        const double t28 = int(t6);
        const double t29 = 0.5 * t13;
        const double t30 = 0.5 * t10 * t17;
        const double t31 = -p0_x + p1_x;
        const double t32 = t6 ? 0 : (2 * t31);
        const double t33 = t10 * t32;
        const double t34 = 0.5 * t19;
        const double t35 = t6 ? (2 * t16) : 0;
        const double t36 = t10 * t35;
        const double t37 = 2 * t7;
        const double t38 = t18 * t37;
        const double t39 = t4 * t8;
        const double t40 = t0 * t13;
        const double t41 = t39 - t40;
        const double t42 = t25 ? 0 : t7;
        const double t43 = t42 * t7;
        const double t44 = t25 ? t16 : t0;
        const double t45 = t31 * t44;
        const double t46 = -t43 + t45;
        const double t47 = t13 * t2;
        const double t48 = t17 * t4;
        const double t49 = t47 - t48;
        const double t50 = t41 * t41 + t46 * t46 + t49 * t49;
        const double t51 = 1.0 / std::sqrt(t50);
        const double t52 = t15 * t4;
        const double t53 = 1.0 / t50;
        const double t54 = t25 ? t2 : 0;
        const double t55 = -t16 * t42 + t31 * t54;
        const double t56 = t54 * t55;
        const double t57 = int(!t25);
        const double t58 = t16 * t44 - t54 * t7;
        const double t59 = t16 * t58;
        const double t60 = t43 - t45;
        const double t61 = t53 * (-t56 + t57 * t59 + t60 * (t0 * t57 + t44));
        const double t62 = -t39 + t40;
        const double t63 = -t0 * t17 + t2 * t8;
        const double t64 = -t47 + t48;
        const double t65 = t62 * t62 + t63 * t63 + t64 * t64;
        const double t66 = 1.0 / std::sqrt(t65);
        const double t67 = t13 * t41;
        const double t68 = t0 * t15 + t17;
        const double t69 = 1.0 / t65;
        const double t70 = t63 * t69;
        const double t71 = t42 * t55;
        const double t72 = -int(t25);
        const double t73 = t0 * t60;
        const double t74 = -t44;
        const double t75 = t53 * (t58 * (t16 * t72 + t74) + t71 + t72 * t73);
        const double t76 = t17 + t22 * t4;
        const double t77 = t0 * t22;
        const double t78 = t41 * t8;
        const double t79 = int(t25);
        const double t80 = -int(!t25);
        const double t81 = t49 * t53;
        const double t82 = t13 + t2 * t28;
        const double t83 = t2 * t23;
        const double t84 = t23 * t4;
        const double t85 = t0 * t28;
        const double t86 = t84 - t85;
        const double t87 = t41 * t86 + t46 * (t8 + t83) + t49 * t82;
        const double t88 = t62 * t69;
        const double t89 = t56 + t59 * t80 + t60 * (t0 * t80 + t74);
        const double t90 = t0 * t23 - t17;
        const double t91 = t41 * t53;
        const double t92 = t17 - t28 * t4;
        const double t93 = -t46 * t85 + t49 * t92 - t78;
        const double t94 = -t13 + t2 * t22;
        const double t95 = t52 - t77;
        const double t96 = -t15 * t2 + t8;
        const double t97 = t41 * t95 - t46 * t96 + t49 * t94;
        J[0] = -t12 * t8;
        J[1] = -t12 * t13;
        J[2] = t14 * (-t11 * t19 + t15);
        J[3] = -t21 * t8;
        J[4] = -t13 * t21;
        J[5] = t14 * (-t19 * t20 + t22);
        J[6] = t14 * (t23 - t24 * t27);
        J[7] = t14 * (-t27 * t29 + t28);
        J[8] = -t30 * t26;
        J[9] = -t24 * t33;
        J[10] = -t29 * t33;
        J[11] = t14 * (t23 - t32 * t34);
        J[12] = -t24 * t36;
        J[13] = -t29 * t36;
        J[14] = t14 * (t28 - t34 * t35);
        J[15] = t14 * (t15 - t24 * t38);
        J[16] = t14 * (t22 - t29 * t38);
        J[17] = -t30 * t37;
        J[18] = -t51 * (t49 * t61 + t52);
        J[19] = t66 * (t68 - t70 * (t46 * t68 + t49 * t52 + t67));
        J[20] = -t51 * (t13 + t41 * t61);
        J[21] = -t51 * (t49 * t75 + t76);
        J[22] = t66 * (t70 * (-t46 * t77 - t49 * t76 + t78) + t77);
        J[23] = t51 * (-t41 * t75 + t8);
        J[24] = t51
            * (-t81
                   * (t55 * (t31 * t79 + t4 * t80) + t58 * (t2 * t79 + t54)
                      + t60 * (-t42 + t7 * t80))
               + t82);
        J[25] = t66 * (t70 * t87 - t8 - t83);
        J[26] = t66 * (t86 + t87 * t88);
        J[27] = -t51 * (t81 * t89 + t84);
        J[28] = t66 * (t70 * (-t46 * t90 - t49 * t84 + t67) + t90);
        J[29] = t51 * (t13 - t89 * t91);
        J[30] = t66 * (t64 * t69 * t93 + t92);
        J[31] = t66 * (t70 * t93 + t85);
        J[32] = -t51 * (t8 + t91 * (t58 * (t16 * t79 + t44) - t71 + t73 * t79));
        J[33] = t51
            * (-t81
                   * (t55 * (t31 * t72 + t4 * t57) + t58 * (t2 * t72 - t54)
                      + t60 * (t42 + t57 * t7))
               + t94);
        J[34] = t66 * (t70 * t97 + t96);
        J[35] = t66 * (t88 * t97 + t95);
    }

    // J is (12×1) flattened in column-major order
    void point_edge_tangent_basis_2D_jacobian(
        double p_x,
        double p_y,
        double e0_x,
        double e0_y,
        double e1_x,
        double e1_y,
        double J[12])
    {
        const double t0 = e0_x - e1_x;
        const double t1 = t0 * t0;
        const double t2 = e0_y - e1_y;
        const double t3 = t2 * t2;
        const double t4 = t1 + t3;
        const double t5 = 1.0 / std::sqrt(t4);
        const double t6 = 1.0 / t4;
        const double t7 = t1 * t6;
        const double t8 = t0 * t2 / pow_1_5(t4);
        const double t9 = t3 * t6;
        const double t10 = -t8;
        J[0] = 0;
        J[1] = 0;
        J[2] = 0;
        J[3] = 0;
        J[4] = t5 * (t7 - 1);
        J[5] = t8;
        J[6] = t8;
        J[7] = t5 * (t9 - 1);
        J[8] = t5 * (1 - t7);
        J[9] = t10;
        J[10] = t10;
        J[11] = t5 * (1 - t9);
    }

    // J is (27×2) flattened in column-major order
    void point_edge_tangent_basis_3D_jacobian(
        double p_x,
        double p_y,
        double p_z,
        double e0_x,
        double e0_y,
        double e0_z,
        double e1_x,
        double e1_y,
        double e1_z,
        double J[54])
    {
        const double t0 = -e1_x;
        const double t1 = e0_x + t0;
        const double t2 = t1 * t1;
        const double t3 = -e1_y;
        const double t4 = e0_y + t3;
        const double t5 = t4 * t4;
        const double t6 = -e1_z;
        const double t7 = e0_z + t6;
        const double t8 = t7 * t7;
        const double t9 = t2 + t5 + t8;
        const double t10 = 1.0 / std::sqrt(t9);
        const double t11 = 1.0 / t9;
        const double t12 = t11 * t2;
        const double t13 = 1.0 / pow_1_5(t9);
        const double t14 = t1 * t13;
        const double t15 = t14 * t4;
        const double t16 = t14 * t7;
        const double t17 = t11 * t5;
        const double t18 = t13 * t4 * t7;
        const double t19 = t11 * t8;
        const double t20 = -t15;
        const double t21 = -t16;
        const double t22 = -t18;
        const double t23 = -p_y;
        const double t24 = e0_y + t23;
        const double t25 = -p_x;
        const double t26 = e0_x + t25;
        const double t27 = t1 * t24 - t26 * t4;
        const double t28 = -e0_x;
        const double t29 = e1_x + t28;
        const double t30 = -e0_z;
        const double t31 = p_z + t30;
        const double t32 = t29 * t31;
        const double t33 = p_x + t28;
        const double t34 = e1_z + t30;
        const double t35 = t33 * t34;
        const double t36 = t32 - t35;
        const double t37 = t27 * t4 + t36 * t7;
        const double t38 = -p_z;
        const double t39 = e0_z + t38;
        const double t40 = -t24 * t7 + t39 * t4;
        const double t41 = t1 * t39 - t26 * t7;
        const double t42 = t27 * t27 + t40 * t40;
        const double t43 = t41 * t41 + t42;
        const double t44 = 1.0 / pow_1_5(t43);
        const double t45 = t40 * t44;
        const double t46 = t36 * t36 + t42;
        const double t47 = 1 / std::sqrt(t46);
        const double t48 = -e0_y;
        const double t49 = p_y + t48;
        const double t50 = e1_y + t48;
        const double t51 = t29 * t49 - t33 * t50;
        const double t52 = -t32 + t35;
        const double t53 = 1.0 / t46;
        const double t54 = t36 * t53;
        const double t55 = 1.0 / std::sqrt(t43);
        const double t56 = 1.0 / t43;
        const double t57 = t27 * t56;
        const double t58 = -t1 * t27 + t40 * t7;
        const double t59 = t40 * t56;
        const double t60 = t41 * t44;
        const double t61 = t31 * t50 - t34 * t49;
        const double t62 = t27 * t53;
        const double t63 = t40 * t53;
        const double t64 = t1 * t36 + t4 * t40;
        const double t65 = t41 * t56;
        const double t66 = t27 * t44;
        const double t67 = e1_y + t23;
        const double t68 = e1_z + t38;
        const double t69 = t27 * t67 + t36 * t68;
        const double t70 = p_z + t6;
        const double t71 = e1_x + t25;
        const double t72 = -t27 * t71 + t40 * t68;
        const double t73 = t36 * t71 + t40 * t67;
        const double t74 = t24 * t27 + t36 * t39;
        const double t75 = t26 * t27 - t39 * t40;
        const double t76 = t24 * t40 + t26 * t36;
        J[0] = 0;
        J[1] = 0;
        J[2] = 0;
        J[3] = 0;
        J[4] = 0;
        J[5] = 0;
        J[6] = 0;
        J[7] = 0;
        J[8] = 0;
        J[9] = t10 * (t12 - 1);
        J[10] = t15;
        J[11] = t16;
        J[12] = t15;
        J[13] = t10 * (t17 - 1);
        J[14] = t18;
        J[15] = t16;
        J[16] = t18;
        J[17] = t10 * (t19 - 1);
        J[18] = t10 * (1 - t12);
        J[19] = t20;
        J[20] = t21;
        J[21] = t20;
        J[22] = t10 * (1 - t17);
        J[23] = t22;
        J[24] = t21;
        J[25] = t22;
        J[26] = t10 * (1 - t19);
        J[27] = -t37 * t45;
        J[28] = t47 * (t34 + t54 * (t34 * t52 + t4 * t51));
        J[29] = t55 * (-t37 * t57 + t4);
        J[30] = t55 * (-t58 * t59 + t7);
        J[31] = t58 * t60;
        J[32] = -t47 * (t1 + t62 * (t29 * t51 + t61 * t7));
        J[33] = -t47 * (t4 + t63 * (t1 * t52 + t50 * t61));
        J[34] = t55 * (t1 - t64 * t65);
        J[35] = t64 * t66;
        J[36] = -t45 * t69;
        J[37] = t47 * (t54 * (t51 * t67 + t52 * t70) + t70);
        J[38] = t55 * (-t57 * t69 + t67);
        J[39] = t55 * (-t59 * t72 + t68);
        J[40] = t60 * t72;
        J[41] = -t47 * (t62 * (t51 * (p_x + t0) + t61 * t68) + t71);
        J[42] = -t47 * (t63 * (t52 * t71 + t61 * (p_y + t3)) + t67);
        J[43] = t55 * (-t65 * t73 + t71);
        J[44] = t66 * t73;
        J[45] = t45 * t74;
        J[46] = t55 * (t39 - t65 * t74);
        J[47] = -t47 * (t24 + t62 * (t39 * t52 + t49 * t51));
        J[48] = -t47 * (t39 + t63 * (t26 * t51 + t31 * t61));
        J[49] = t60 * t75;
        J[50] = t55 * (t26 - t57 * t75);
        J[51] = t55 * (t24 - t59 * t76);
        J[52] = t47 * (t33 + t54 * (t24 * t61 + t33 * t52));
        J[53] = -t66 * t76;
    }

    // J is (36×2) flattened in column-major order
    void edge_edge_tangent_basis_jacobian(
        double ea0_x,
        double ea0_y,
        double ea0_z,
        double ea1_x,
        double ea1_y,
        double ea1_z,
        double eb0_x,
        double eb0_y,
        double eb0_z,
        double eb1_x,
        double eb1_y,
        double eb1_z,
        double J[72])
    {
        const double t0 = ea0_x - ea1_x;
        const double t1 = t0 * t0;
        const double t2 = ea0_y - ea1_y;
        const double t3 = t2 * t2;
        const double t4 = ea0_z - ea1_z;
        const double t5 = t4 * t4;
        const double t6 = t3 + t5;
        const double t7 = t1 + t6;
        const double t8 = 1 / std::sqrt(t7);
        const double t9 = 1.0 / t7;
        const double t10 = t1 * t9;
        const double t11 = 1.0 / pow_1_5(t7);
        const double t12 = t0 * t2;
        const double t13 = t11 * t12;
        const double t14 = t0 * t4;
        const double t15 = t11 * t14;
        const double t16 = t3 * t9;
        const double t17 = t2 * t4;
        const double t18 = t11 * t17;
        const double t19 = t5 * t9;
        const double t20 = -t13;
        const double t21 = -t15;
        const double t22 = -t18;
        const double t23 = -ea0_x + ea1_x;
        const double t24 = -eb0_z + eb1_z;
        const double t25 = t23 * t24;
        const double t26 = -ea0_z + ea1_z;
        const double t27 = -eb0_x + eb1_x;
        const double t28 = t26 * t27;
        const double t29 = t25 - t28;
        const double t30 = eb0_z - eb1_z;
        const double t31 = t2 * t30;
        const double t32 = eb0_y - eb1_y;
        const double t33 = t32 * t4;
        const double t34 = -t33;
        const double t35 = t31 + t34;
        const double t36 = t0 * t29 + t2 * t35;
        const double t37 = t0 * t32;
        const double t38 = eb0_x - eb1_x;
        const double t39 = t2 * t38;
        const double t40 = -t39;
        const double t41 = t37 + t40;
        const double t42 = t2 * t41 + t29 * t4;
        const double t43 = t0 * t41;
        const double t44 = t35 * t4;
        const double t45 = t43 - t44;
        const double t46 = t36 * t36 + t42 * t42 + t45 * t45;
        const double t47 = 1 / std::sqrt(t46);
        const double t48 = 1.0 / t46;
        const double t49 = 2 * t37;
        const double t50 = t39 - t49;
        const double t51 = -t43 + t44;
        const double t52 = t2 * t32;
        const double t53 = t30 * t4;
        const double t54 = t52 + t53;
        const double t55 = -ea0_y + ea1_y;
        const double t56 = -eb0_y + eb1_y;
        const double t57 = t23 * t56;
        const double t58 = t27 * t55;
        const double t59 = t57 - t58;
        const double t60 = t55 * t59;
        const double t61 = -t25;
        const double t62 = t28 + t61;
        const double t63 = t26 * t62;
        const double t64 = t60 - t63;
        const double t65 = t54 * t64;
        const double t66 = t38 * t4;
        const double t67 = t0 * t30;
        const double t68 = 2 * t67;
        const double t69 = t66 - t68;
        const double t70 = t23 * t62;
        const double t71 = t24 * t55;
        const double t72 = t26 * t56;
        const double t73 = -t72;
        const double t74 = t71 + t73;
        const double t75 = t55 * t74;
        const double t76 = t70 - t75;
        const double t77 = t48 * (t50 * t51 - t65 - t69 * t76);
        const double t78 = t0 * t38;
        const double t79 = t53 + t78;
        const double t80 = t51 * t79;
        const double t81 = 2 * t31;
        const double t82 = t33 - t81;
        const double t83 = 2 * t39;
        const double t84 = t37 - t83;
        const double t85 = t48 * (-t64 * t84 - t76 * t82 + t80);
        const double t86 = t52 + t78;
        const double t87 = t76 * t86;
        const double t88 = 2 * t33;
        const double t89 = t31 - t88;
        const double t90 = -2 * t66 + t67;
        const double t91 = t48 * (t51 * t89 - t64 * t90 - t87);
        const double t92 = t51 * t51 + t64 * t64 + t76 * t76;
        const double t93 = 1.0 / std::sqrt(t92);
        const double t94 = -t70 + t75;
        const double t95 = t23 * t59 - t26 * t74;
        const double t96 = -t60 + t63;
        const double t97 = 1.0 / t92;
        const double t98 = t51 * t97;
        const double t99 = t40 + t49;
        const double t100 = -t66 + t68;
        const double t101 = t48 * (-t100 * t76 + t51 * t99 + t65);
        const double t102 = t34 + t81;
        const double t103 = -t37 + t83;
        const double t104 = t48 * (t102 * t76 + t103 * t64 + t80);
        const double t105 = 2 * t28;
        const double t106 = t105 + t61;
        const double t107 =
            t106 * t96 + t94 * (t0 * t27 + t32 * t55) + t95 * (t33 - t71 + t72);
        const double t108 = t64 * t97;
        const double t109 = -t31 + t88;
        const double t110 = t36 * t48;
        const double t111 = t12 * t51;
        const double t112 = t14 * t76;
        const double t113 = t6 * t64;
        const double t114 = t111 - t112 + t113;
        const double t115 = t114 * t48;
        const double t116 = t1 + t5;
        const double t117 = t116 * t51 + t12 * t64 + t17 * t76;
        const double t118 = t42 * t48;
        const double t119 = t26 * t94;
        const double t120 = t17 * t51;
        const double t121 = t14 * t64;
        const double t122 = t1 + t3;
        const double t123 = t122 * t76;
        const double t124 = t120 - t121 + t123;
        const double t125 = t2 * t23;
        const double t126 = t0 * t119 + t125 * t95 + t96 * (t26 * t26 + t3);
        const double t127 = -t14;
        const double t128 = t76 * t97;
        const double t129 = t4 * t55;
        const double t130 = t125 * t96 + t129 * t94 + t95 * (t23 * t23 + t5);
        const double t131 =
            t0 * t26 * t96 + t129 * t95 + t94 * (t1 + t55 * t55);
        J[0] = t8 * (t10 - 1);
        J[1] = t13;
        J[2] = t15;
        J[3] = t13;
        J[4] = t8 * (t16 - 1);
        J[5] = t18;
        J[6] = t15;
        J[7] = t18;
        J[8] = t8 * (t19 - 1);
        J[9] = t8 * (1 - t10);
        J[10] = t20;
        J[11] = t21;
        J[12] = t20;
        J[13] = t8 * (1 - t16);
        J[14] = t22;
        J[15] = t21;
        J[16] = t22;
        J[17] = t8 * (1 - t19);
        J[18] = 0;
        J[19] = 0;
        J[20] = 0;
        J[21] = 0;
        J[22] = 0;
        J[23] = 0;
        J[24] = 0;
        J[25] = 0;
        J[26] = 0;
        J[27] = 0;
        J[28] = 0;
        J[29] = 0;
        J[30] = 0;
        J[31] = 0;
        J[32] = 0;
        J[33] = 0;
        J[34] = 0;
        J[35] = 0;
        J[36] = t47 * (-t42 * t77 + t54);
        J[37] = t47 * (t45 * t77 + t50);
        J[38] = t47 * (t36 * t77 + t69);
        J[39] = t47 * (-t42 * t85 + t84);
        J[40] = t47 * (t45 * t85 + t79);
        J[41] = t47 * (t36 * t85 + t82);
        J[42] = t47 * (-t42 * t91 + t90);
        J[43] = t93
            * (t89
               - t98
                   * (t94 * (t55 * t56 + t78) + t95 * (t4 * t56 + t74)
                      + t96 * (t26 * t38 + t29)));
        J[44] = t47 * (t36 * t91 + t86);
        J[45] = -t47 * (t101 * t42 + t54);
        J[46] = t47 * (t101 * t45 + t99);
        J[47] = t47 * (t100 + t101 * t36);
        J[48] = t47 * (t103 + t104 * t42);
        J[49] = -t93
            * (t79
               + t98
                   * (t94 * (2 * t71 + t73) + t95 * (t23 * t38 + t24 * t4)
                      + t96 * (t39 - t57 + t58)));
        J[50] = t47 * (t102 - t104 * t36);
        J[51] = t93 * (t106 + t107 * t108);
        J[52] = t93 * (-t107 * t98 + t109);
        J[53] =
            t47 * (t110 * (t109 * t51 + t64 * (-t105 + t25) + t87) - t52 - t78);
        J[54] = -t47 * (t115 * t42 + t6);
        J[55] = t47 * (t115 * t45 + t12);
        J[56] = t47 * (t110 * t114 + t14);
        J[57] = t47 * (t117 * t118 + t12);
        J[58] = -t93
            * (t116
               + t98 * (t119 * t55 + t12 * t96 + t95 * (t0 * t23 + t26 * t4)));
        J[59] = t47 * (-t110 * t117 + t17);
        J[60] = t47 * (-t118 * t124 + t14);
        J[61] = t47 * (t124 * t45 * t48 + t17);
        J[62] = t47 * (-t1 + t110 * t124 - t3);
        J[63] = t47 * (-t118 * (-t111 + t112 - t113) + t6);
        J[64] = -t93 * (t12 + t126 * t98);
        J[65] = t93 * (t126 * t128 + t127);
        J[66] = t93 * (t108 * t130 - t12);
        J[67] = t93 * (t116 - t130 * t98);
        J[68] = t93 * (t128 * t130 - t17);
        J[69] = t93 * (t108 * t131 + t127);
        J[70] = -t93 * (t131 * t98 + t17);
        J[71] = t47 * (t110 * (-t120 + t121 - t123) + t122);
    }

    // J is (36×2) flattened in column-major order
    void point_triangle_tangent_basis_jacobian(
        double p_x,
        double p_y,
        double p_z,
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double J[72])
    {
        const double t0 = t0_x - t1_x;
        const double t1 = t0 * t0;
        const double t2 = -t1_y;
        const double t3 = t0_y + t2;
        const double t4 = t3 * t3;
        const double t5 = -t1_z;
        const double t6 = t0_z + t5;
        const double t7 = t6 * t6;
        const double t8 = t4 + t7;
        const double t9 = t1 + t8;
        const double t10 = 1.0 / std::sqrt(t9);
        const double t11 = 1.0 / t9;
        const double t12 = t1 * t11;
        const double t13 = 1.0 / pow_1_5(t9);
        const double t14 = t0 * t3;
        const double t15 = t13 * t14;
        const double t16 = t0 * t6;
        const double t17 = t13 * t16;
        const double t18 = t11 * t4;
        const double t19 = t3 * t6;
        const double t20 = t13 * t19;
        const double t21 = t11 * t7;
        const double t22 = -t15;
        const double t23 = -t17;
        const double t24 = -t20;
        const double t25 = -t0_x;
        const double t26 = t1_x + t25;
        const double t27 = -t0_z;
        const double t28 = t27 + t2_z;
        const double t29 = t26 * t28;
        const double t30 = t25 + t2_x;
        const double t31 = t1_z + t27;
        const double t32 = t30 * t31;
        const double t33 = t29 - t32;
        const double t34 = -t2_z;
        const double t35 = t0_z + t34;
        const double t36 = t3 * t35;
        const double t37 = -t2_y;
        const double t38 = t0_y + t37;
        const double t39 = t38 * t6;
        const double t40 = -t39;
        const double t41 = t36 + t40;
        const double t42 = t0 * t33 + t3 * t41;
        const double t43 = t0 * t38;
        const double t44 = -t2_x;
        const double t45 = t0_x + t44;
        const double t46 = t3 * t45;
        const double t47 = -t46;
        const double t48 = t43 + t47;
        const double t49 = t3 * t48 + t33 * t6;
        const double t50 = t0 * t48;
        const double t51 = t41 * t6;
        const double t52 = t50 - t51;
        const double t53 = t42 * t42 + t49 * t49 + t52 * t52;
        const double t54 = 1.0 / std::sqrt(t53);
        const double t55 = 1.0 / t53;
        const double t56 = t1_y + t37;
        const double t57 = t3 * t56;
        const double t58 = t1_z + t34;
        const double t59 = t58 * t6;
        const double t60 = t57 + t59;
        const double t61 = -t0_y;
        const double t62 = t1_y + t61;
        const double t63 = t2_y + t61;
        const double t64 = t26 * t63;
        const double t65 = t30 * t62;
        const double t66 = t64 - t65;
        const double t67 = t62 * t66;
        const double t68 = -t29;
        const double t69 = t32 + t68;
        const double t70 = t31 * t69;
        const double t71 = t67 - t70;
        const double t72 = -t43;
        const double t73 = t46 + t72;
        const double t74 = -t0 * t56 + t73;
        const double t75 = -t50 + t51;
        const double t76 = t45 * t6;
        const double t77 = t0 * t35;
        const double t78 = t76 - t77;
        const double t79 = -t0 * t58 + t78;
        const double t80 = t26 * t69;
        const double t81 = t28 * t62;
        const double t82 = t31 * t63;
        const double t83 = -t82;
        const double t84 = t81 + t83;
        const double t85 = t62 * t84;
        const double t86 = t80 - t85;
        const double t87 = t55 * (-t60 * t71 + t74 * t75 - t79 * t86);
        const double t88 = t71 * t71 + t75 * t75 + t86 * t86;
        const double t89 = 1.0 / std::sqrt(t88);
        const double t90 = t2_z + t5;
        const double t91 = -t67 + t70;
        const double t92 = -t64 + t65;
        const double t93 = t26 * t66 - t31 * t84;
        const double t94 = -t80 + t85;
        const double t95 = 1.0 / t88;
        const double t96 = t75 * t95;
        const double t97 = t1_x + t44;
        const double t98 = t3 * t97;
        const double t99 = t0 * t97;
        const double t100 = t59 + t99;
        const double t101 = -t36;
        const double t102 = t101 - t3 * t58 + t39;
        const double t103 = t55 * (t100 * t75 - t102 * t86 - t71 * (t48 - t98));
        const double t104 = t6 * t97;
        const double t105 = t57 + t99;
        const double t106 = t41 - t56 * t6;
        const double t107 = -t76;
        const double t108 =
            t55 * (-t105 * t86 + t106 * t75 - t71 * (-t104 + t107 + t77));
        const double t109 = t2 + t2_y;
        const double t110 = t3 * t38;
        const double t111 = t35 * t6;
        const double t112 = t110 + t111;
        const double t113 = 2 * t43 + t47;
        const double t114 = t107 + 2 * t77;
        const double t115 = t55 * (t112 * t71 + t113 * t75 - t114 * t86);
        const double t116 = t0 * t45;
        const double t117 = t111 + t116;
        const double t118 = 2 * t36 + t40;
        const double t119 = 2 * t46 + t72;
        const double t120 = t55 * (t117 * t75 + t118 * t86 + t119 * t71);
        const double t121 = 2 * t32;
        const double t122 = t121 + t68;
        const double t123 =
            t122 * t91 + t93 * (t39 - t81 + t82) + t94 * (t0 * t30 + t38 * t62);
        const double t124 = t71 * t95;
        const double t125 = t101 + 2 * t39;
        const double t126 = t42 * t55;
        const double t127 = t26 * t3;
        const double t128 = t0 * t31;
        const double t129 = t127 * t93 + t128 * t94 + t91 * (t31 * t31 + t4);
        const double t130 = -t16;
        const double t131 = t86 * t95;
        const double t132 = t6 * t62;
        const double t133 = t127 * t91 + t132 * t94 + t93 * (t26 * t26 + t7);
        const double t134 = t128 * t91 + t132 * t93 + t94 * (t1 + t62 * t62);
        const double t135 = t1 + t4;
        J[0] = 0;
        J[1] = 0;
        J[2] = 0;
        J[3] = 0;
        J[4] = 0;
        J[5] = 0;
        J[6] = 0;
        J[7] = 0;
        J[8] = 0;
        J[9] = t10 * (t12 - 1);
        J[10] = t15;
        J[11] = t17;
        J[12] = t15;
        J[13] = t10 * (t18 - 1);
        J[14] = t20;
        J[15] = t17;
        J[16] = t20;
        J[17] = t10 * (t21 - 1);
        J[18] = t10 * (1 - t12);
        J[19] = t22;
        J[20] = t23;
        J[21] = t22;
        J[22] = t10 * (1 - t18);
        J[23] = t24;
        J[24] = t23;
        J[25] = t24;
        J[26] = t10 * (1 - t21);
        J[27] = 0;
        J[28] = 0;
        J[29] = 0;
        J[30] = 0;
        J[31] = 0;
        J[32] = 0;
        J[33] = 0;
        J[34] = 0;
        J[35] = 0;
        J[36] = 0;
        J[37] = 0;
        J[38] = 0;
        J[39] = 0;
        J[40] = 0;
        J[41] = 0;
        J[42] = 0;
        J[43] = 0;
        J[44] = 0;
        J[45] = t54 * (-t49 * t87 + t60);
        J[46] = t89
            * (t74
               - t96
                   * (t91 * (t31 * t90 + t57) + t93 * (t26 * t56 + t92)
                      + t94 * (t0 * t90 + t69)));
        J[47] = t54 * (t42 * t87 + t79);
        J[48] = -t54 * (t103 * t49 + t73 + t98);
        J[49] = t54 * (t100 + t103 * t52);
        J[50] = t54 * (t102 + t103 * t42);
        J[51] = -t54 * (t104 + t108 * t49 + t78);
        J[52] = t89
            * (t106
               - t96
                   * (t91 * (t31 * t97 + t33) + t93 * (t109 * t6 + t84)
                      + t94 * (t109 * t62 + t99)));
        J[53] = t54 * (t105 + t108 * t42);
        J[54] = -t54 * (t112 + t115 * t49);
        J[55] = t54 * (t113 + t115 * t52);
        J[56] = t54 * (t114 + t115 * t42);
        J[57] = t54 * (t119 + t120 * t49);
        J[58] = -t89
            * (t117
               + t96
                   * (t91 * (t46 + t92) + t93 * (t26 * t45 + t28 * t6)
                      + t94 * (2 * t81 + t83)));
        J[59] = t54 * (t118 - t120 * t42);
        J[60] = t89 * (t122 + t123 * t124);
        J[61] = t89 * (-t123 * t96 + t125);
        J[62] = t54
            * (-t110 - t116
               + t126
                   * (t125 * t75 + t71 * (-t121 + t29) + t86 * (t110 + t116)));
        J[63] = t54 * (-t49 * t55 * (-t14 * t75 + t16 * t86 - t71 * t8) + t8);
        J[64] = -t89 * (t129 * t96 + t14);
        J[65] = t89 * (t129 * t131 + t130);
        J[66] = t89 * (t124 * t133 - t14);
        J[67] = t89 * (t1 - t133 * t96 + t7);
        J[68] = t89 * (t131 * t133 - t19);
        J[69] = t89 * (t124 * t134 + t130);
        J[70] = -t89 * (t134 * t96 + t19);
        J[71] = t54 * (t126 * (-t135 * t86 + t16 * t71 - t19 * t75) + t135);
    }

} // namespace autogen
} // namespace ipc
