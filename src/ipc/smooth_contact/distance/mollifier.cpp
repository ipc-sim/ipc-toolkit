#include "mollifier.hpp"
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>

namespace ipc {

/// @brief Compute the gradient of the mollifier function wrt. 4 edge points and the distance squared
std::pair<double, Vector<double, 13>> edge_edge_mollifier_grad(
    const Vector3<double> &ea0, const Vector3<double> &ea1,
    const Vector3<double> &eb0, const Vector3<double> &eb1,
    const double &dist_sqr)
{
    const double da0 = point_edge_distance(ea0, eb0, eb1);
    const double da1 = point_edge_distance(ea1, eb0, eb1);
    const double db0 = point_edge_distance(eb0, ea0, ea1);
    const double db1 = point_edge_distance(eb1, ea0, ea1);

    Vector12d da0_wrt_x, da1_wrt_x, db0_wrt_x, db1_wrt_x;
    {
        Vector9d tmp = point_edge_distance_gradient(ea0, eb0, eb1);
        da0_wrt_x << tmp.head(3), Eigen::Vector3d::Zero(), tmp.tail(6);

        tmp = point_edge_distance_gradient(ea1, eb0, eb1);
        da1_wrt_x << Eigen::Vector3d::Zero(), tmp;

        tmp = point_edge_distance_gradient(eb0, ea0, ea1);
        db0_wrt_x << tmp.tail(6), tmp.head(3), Eigen::Vector3d::Zero();

        tmp = point_edge_distance_gradient(eb1, ea0, ea1);
        db1_wrt_x << tmp.tail(6), Eigen::Vector3d::Zero(), tmp.head(3);
    }

    const double la = point_point_distance(ea0, ea1);
    const double lb = point_point_distance(eb0, eb1);

    const double db = lb * mollifier_threshold_eps;
    const double da = la * mollifier_threshold_eps;

    Vector6d da_wrt_x = point_point_distance_gradient(ea0, ea1) * mollifier_threshold_eps;
    Vector6d db_wrt_x = point_point_distance_gradient(eb0, eb1) * mollifier_threshold_eps;

    const double ma0 = (da0 - dist_sqr) / db;
    const double ma1 = (da1 - dist_sqr) / db;
    const double mb0 = (db0 - dist_sqr) / da;
    const double mb1 = (db1 - dist_sqr) / da;

    const double a = Math<double>::mollifier(ma0);
    const double b = Math<double>::mollifier(ma1);
    const double c = Math<double>::mollifier(mb0);
    const double d = Math<double>::mollifier(mb1);

    const double a_wrt_ma0 = Math<double>::mollifier_grad(ma0);
    const double b_wrt_ma1 = Math<double>::mollifier_grad(ma1);
    const double c_wrt_mb0 = Math<double>::mollifier_grad(mb0);
    const double d_wrt_mb1 = Math<double>::mollifier_grad(mb1);

    const double ab = a * b, cd = c * d;
    Vector12d ab_grad = a * b_wrt_ma1 * da1_wrt_x + b * a_wrt_ma0 * da0_wrt_x;
    ab_grad.tail(6)  -= db_wrt_x * (a * b_wrt_ma1 * ma1 + b * a_wrt_ma0 * ma0);
    ab_grad /= db;
    Vector12d cd_grad = c * d_wrt_mb1 * db1_wrt_x + d * c_wrt_mb0 * db0_wrt_x;
    cd_grad.head(6)  -= da_wrt_x * (c * d_wrt_mb1 * mb1 + d * c_wrt_mb0 * mb0);
    cd_grad /= da;

    Vector<double, 13> out;
    out.head(12) = ab_grad * cd + cd_grad * ab;
    out(12) = - (a * b_wrt_ma1 + b * a_wrt_ma0) * cd / db - (c * d_wrt_mb1 + d * c_wrt_mb0) * ab / da;

    return std::make_pair(ab * cd, out);
}

/// @brief Compute the hessian of the mollifier function wrt. 4 edge points and the distance squared
Eigen::Matrix<double, 13, 13> edge_edge_mollifier_hessian(
    const Vector3<double> &ea0, const Vector3<double> &ea1,
    const Vector3<double> &eb0, const Vector3<double> &eb1, 
    const std::array<HEAVISIDE_TYPE, 4> mtypes,
    const double &dist_sqr)
{
    return Eigen::Matrix<double, 13, 13>::Zero();
}

}
