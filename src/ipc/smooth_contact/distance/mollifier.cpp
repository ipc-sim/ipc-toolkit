#include "mollifier.hpp"

namespace ipc {

/// @brief Compute the gradient of the mollifier function wrt. 4 edge points and the distance squared
Vector<double, 13> edge_edge_mollifier_grad(
    const Vector3<double> &ea0, const Vector3<double> &ea1,
    const Vector3<double> &eb0, const Vector3<double> &eb1, 
    const std::array<HEAVISIDE_TYPE, 4> mtypes,
    const double &dist_sqr)
{
    return Vector<double, 13>::Zero();
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
