#pragma once

#include <ipc/config.hpp>
#include "eigen_ext.hpp"

namespace ipc {
    template <typename scalar>
    scalar abs(const scalar &x);

    template <typename scalar>
    scalar sqr(const scalar &x)
    {
        return x * x;
    }

    template <typename scalar>
    scalar cubic(const scalar &x)
    {
        return x * x * x;
    }

    template <typename scalar>
    scalar cubic_spline(const scalar &x);

    /// @brief support is [-1, 1]
    /// @tparam scalar 
    /// @param x 
    /// @return 
    template <typename scalar>
    scalar quadratic_spline(const scalar &x);

    enum class HEAVISIDE_TYPE {
        ZERO = 0, ONE = 1, VARIANT = 2
    };

    struct ORIENTATION_TYPES {

        static inline HEAVISIDE_TYPE compute_type(const double &val, const double &alpha, const double &beta);

        int size() const { return size_; }
        inline void set_size(const int size);
        HEAVISIDE_TYPE tangent_type(const int &i) const { return tangent_types[i]; }
        HEAVISIDE_TYPE normal_type(const int &i) const { return normal_types[i]; }
        HEAVISIDE_TYPE& tangent_type(const int &i) { return tangent_types[i]; }
        HEAVISIDE_TYPE& normal_type(const int &i) { return normal_types[i]; }

        int size_ = 0;
        std::vector<HEAVISIDE_TYPE> tangent_types, normal_types;
    };

    template <typename scalar>
    scalar smooth_heaviside(const scalar &x, const double alpha, const double beta = 0);

    constexpr double mollifier_threshold_eps = 1e-3;

    template <typename scalar>
    scalar mollifier(const scalar &x);

    // support is [0, 1]
    template <typename scalar>
    scalar inv_barrier(const scalar &x, const double &r);

    template <typename scalar>
    scalar L_ns(const scalar &x);

    template <typename scalar>
    scalar cross2(const Eigen::Ref<const Vector2<scalar>> &a, const Eigen::Ref<const Vector2<scalar>> &b);

    // linear solve for 2x2 matrix
    template <typename scalar>
    Vector2<scalar> linear_solve(const Eigen::Ref<const Matrix2<scalar>> &A, const Eigen::Ref<const Vector2<scalar>> &b);

    template <class T, int rows, int dim, int max_rows=rows>
    Eigen::Matrix<T, rows, dim, Eigen::ColMajor, max_rows, dim> slice_positions(const Eigen::VectorXd &positions);

    enum class FD_RULE { CENTRAL, LEFT, RIGHT };
    
    void my_finite_gradient(const Eigen::VectorXd& x, const std::function<double(const Eigen::VectorXd&)> &f, Eigen::VectorXd &grad, FD_RULE rule = FD_RULE::CENTRAL, const double eps = 1e-7);
}

#include "math.tpp"
