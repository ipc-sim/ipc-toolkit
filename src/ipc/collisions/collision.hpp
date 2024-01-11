
#pragma once

#include <ipc/candidates/collision_stencil.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

struct ParameterType
{
    ParameterType(const double &_eps, const double &_alpha, const double &_a, const double &_r, const int &_n_quadrature) : 
    eps(_eps), alpha(_alpha), a(_a), r(_r), n_quadrature(_n_quadrature) 
    {
        assert(r > 0);
        assert(a >= 0);
        assert(eps > 0);
        assert(alpha > 0);
        assert(n_quadrature > 0);
    }
    ParameterType() = delete;

    double eps;
    double alpha;
    double a;
    double r;
    int n_quadrature;
};

template <int max_vert = 4>
class Collision : virtual public CollisionStencil<max_vert> {
public:
    Collision() = default;

    Collision(
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient);

    virtual ~Collision() { }

    // -- non distance type potential ----

    virtual double operator()(
        const Vector<double, -1, 3*max_vert>& positions, 
        const ParameterType &params) const { return 0.; }

    virtual Vector<double, -1, 3*max_vert> gradient(
        const Vector<double, -1, 3*max_vert>& positions, 
        const ParameterType &params) const { return Vector<double, -1, 3*max_vert>::Zero(positions.size()); }

    virtual MatrixMax<double, 3*max_vert, 3*max_vert> hessian(
        const Vector<double, -1, 3*max_vert>& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd = false) const { return MatrixMax<double, 3*max_vert, 3*max_vert>::Zero(positions.size(), positions.size()); }

    // -- Distance mollifier ---------------------------------------------------

    /// @brief Does the distance potentially have to be mollified?
    virtual bool is_mollified() const { return false; }

    /// @brief Compute the mollifier threshold for the distance.
    /// @param rest_positions The stencil's rest vertex positions.
    /// @return The mollifier threshold.
    virtual double mollifier_threshold(const Vector<double, -1, 3*max_vert>& rest_positions) const
    {
        return std::numeric_limits<double>::quiet_NaN(); // No mollifier
    }

    /// @brief Compute the mollifier for the distance.
    /// @param positions The stencil's vertex positions.
    /// @return The mollifier value.
    virtual double mollifier(const Vector<double, -1, 3*max_vert>& positions) const;

    /// @brief Compute the mollifier for the distance.
    /// @param positions The stencil's vertex positions.
    /// @param eps_x The mollifier's threshold.
    /// @return The mollifier value.
    virtual double mollifier(const Vector<double, -1, 3*max_vert>& positions, double eps_x) const;

    /// @brief Compute the gradient of the mollifier for the distance wrt the positions.
    /// @param positions The stencil's vertex positions.
    /// @return The mollifier gradient.
    virtual Vector<double, -1, 3*max_vert>
    mollifier_gradient(const Vector<double, -1, 3*max_vert>& positions) const;

    /// @brief Compute the gradient of the mollifier for the distance wrt the positions.
    /// @param positions The stencil's vertex positions.
    /// @param eps_x The mollifier's threshold.
    /// @return The mollifier gradient.
    virtual Vector<double, -1, 3*max_vert>
    mollifier_gradient(const Vector<double, -1, 3*max_vert>& positions, double eps_x) const;

    /// @brief Compute the Hessian of the mollifier for the distance wrt the positions.
    /// @param positions The stencil's vertex positions.
    /// @return The mollifier Hessian.
    virtual MatrixMax<double, 3*max_vert, 3*max_vert> mollifier_hessian(const Vector<double, -1, 3*max_vert>& positions) const;

    /// @brief Compute the Hessian of the mollifier for the distance wrt the positions.
    /// @param positions The stencil's vertex positions.
    /// @param eps_x The mollifier's threshold.
    /// @return The mollifier Hessian.
    virtual MatrixMax<double, 3*max_vert, 3*max_vert>
    mollifier_hessian(const Vector<double, -1, 3*max_vert>& positions, double eps_x) const;

    /// @brief Compute the gradient of the mollifier for the distance w.r.t. rest positions.
    /// @param rest_positions The stencil's rest vertex positions.
    /// @param positions The stencil's vertex positions.
    /// @return The mollifier gradient w.r.t. rest positions.
    virtual Vector12d mollifier_gradient_wrt_x(
        const Vector<double, -1, 3*max_vert>& rest_positions,
        const Vector<double, -1, 3*max_vert>& positions) const;

    /// @brief Compute the jacobian of the distance mollifier's gradient w.r.t. rest positions.
    /// @param rest_positions The stencil's rest vertex positions.
    /// @param positions The stencil's vertex positions.
    /// @return The jacobian of the mollifier's gradient w.r.t. rest positions.
    virtual Matrix12d mollifier_gradient_jacobian_wrt_x(
        const Vector<double, -1, 3*max_vert>& rest_positions,
        const Vector<double, -1, 3*max_vert>& positions) const;

    // -------------------------------------------------------------------------

    /// @brief The minimum separation distance.
    double dmin = 0;

    /// @brief The term's weight (e.g., collision area)
    double weight = 1;

    /// @brief The gradient of the term's weight wrt the rest positions.
    Eigen::SparseVector<double> weight_gradient;
};

} // namespace ipc
