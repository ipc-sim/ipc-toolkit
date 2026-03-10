#pragma once

#include <ipc/math/interval.hpp>

namespace ipc {

/// @brief Compute the sinc function: \f$ \frac{\sin(x)}{x} \f$
/// @param x The value for which to compute the sinc function
/// @return The value of sinc(x)
double sinc(const double x);

#ifdef IPC_TOOLKIT_WITH_FILIB

/// @brief Compute the sinc function: \f$ \frac{\sin(x)}{x} \f$ over an interval
/// @param x The interval for which to compute the sinc function
/// @return The interval of sinc(x)
filib::Interval sinc(const filib::Interval& x);

#endif

/// @brief Compute the sinc of the norm of a vector (\f$\operatorname{sinc}(\|x\|)\f$)
/// @tparam T The type of the elements of the vector
/// @param x The vector for which to compute the sinc of the norm
/// @return The value of \f$\operatorname{sinc}(\|x\|)\f$
template <typename T> T sinc_norm_x(Eigen::ConstRef<Eigen::VectorX<T>> x)
{
    if constexpr (std::is_same_v<T, double>) {
        return sinc(x.norm());
    } else {
        return sinc(norm(x));
    }
}

/// @brief Compute the gradient of the sinc of the norm of a vector
/// @param x The vector for which to compute the gradient of the sinc of the norm
/// @return The gradient of the sinc of the norm of the vector
VectorMax3d sinc_norm_x_grad(Eigen::ConstRef<VectorMax3d> x);

/// @brief Compute the Hessian of the sinc of the norm of a vector
/// @param x The vector for which to compute the Hessian of the sinc of the norm
/// @return The Hessian of the sinc of the norm of the vector
MatrixMax3d sinc_norm_x_hess(Eigen::ConstRef<VectorMax3d> x);

} // namespace ipc