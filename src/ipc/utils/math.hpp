#pragma once

#include <ipc/config.hpp>
#include "eigen_ext.hpp"

namespace ipc {
    enum class HEAVISIDE_TYPE {
        ZERO = 0, ONE = 1, VARIANT = 2
    };

    struct ORIENTATION_TYPES {

        static HEAVISIDE_TYPE compute_type(const double &val, const double &alpha, const double &beta);

        int size() const { return size_; }
        void set_size(const int size);
        HEAVISIDE_TYPE tangent_type(const int &i) const { return tangent_types[i]; }
        HEAVISIDE_TYPE normal_type(const int &i) const { return normal_types[i]; }
        HEAVISIDE_TYPE& tangent_type(const int &i) { return tangent_types[i]; }
        HEAVISIDE_TYPE& normal_type(const int &i) { return normal_types[i]; }

        int size_ = 0;
        std::vector<HEAVISIDE_TYPE> tangent_types, normal_types;
    };

    constexpr double mollifier_threshold_eps = 1e-3;

    template <typename scalar>
    struct Math {
        static double sign(const double &x);
        static scalar abs(const scalar &x);
        static scalar sqr(const scalar &x);
        static scalar cubic(const scalar &x);

        static scalar cubic_spline(const scalar &x);
        static double cubic_spline_grad(const double &x);
        static double cubic_spline_hess(const double &x);

        /// @brief support is [-1, 1]
        /// @tparam scalar 
        /// @param x 
        /// @return 
        static scalar quadratic_spline(const scalar &x);

        static scalar smooth_heaviside(const scalar &x, const double alpha, const double beta = 0);
        static double smooth_heaviside_grad(const double &x, const double alpha, const double beta = 0);
        static double smooth_heaviside_hess(const double &x, const double alpha, const double beta = 0);

        static scalar mollifier(const scalar &x);
        static double mollifier_grad(const double &x);
        static double mollifier_hess(const double &x);

        // support is [0, 1]
        static scalar inv_barrier(const scalar &x, const double &r);
        static double inv_barrier_grad(const double &x, const double &r);
        static double inv_barrier_hess(const double &x, const double &r);

        static scalar L_ns(const scalar &x);

        static scalar cross2(const Eigen::Ref<const Vector2<scalar>> &a, const Eigen::Ref<const Vector2<scalar>> &b);
    };

    // gradient is symmetric
    std::tuple<Eigen::Vector3d, Eigen::Matrix3d> normalize_vector_grad(const Eigen::Ref<const Eigen::Vector3d> &t);
    // hessian is symmetric wrt. the three dimensions
    std::tuple<Eigen::Vector3d, Eigen::Matrix3d, std::array<Eigen::Matrix<double, 3, 3>, 3>> normalize_vector_hess(const Eigen::Ref<const Eigen::Vector3d> &t);
    
    template <class T, int rows, int cols, int max_rows=rows>
    inline Eigen::Matrix<T, rows, cols, (max_rows > 1 ? Eigen::ColMajor : Eigen::RowMajor), max_rows, cols> slice_positions(const Eigen::VectorXd &positions, const int offset = 0)
    {
        assert(cols > 0);
        const int nrows = rows > 0 ? rows : positions.size() / cols;
        Eigen::Matrix<T, rows, cols, (max_rows > 1 ? Eigen::ColMajor : Eigen::RowMajor), max_rows, cols> points;
        points.setZero(nrows, cols);
        
        for (int i = 0, id = 0; i < nrows; i++)
            for (int d = 0; d < cols; d++, id++)
                if constexpr (std::is_same<T, double>::value)
                    points(i, d) = positions(id);
                else
                    points(i, d) = T(id + offset, positions(id));

        return points;
    }

    template <class T, int rows, int cols, int max_rows=rows>
    inline Eigen::Matrix<double, rows, cols, (max_rows > 1 ? Eigen::ColMajor : Eigen::RowMajor), max_rows, cols> autodiff_to_double(const Eigen::Matrix<T, rows, cols, (max_rows > 1 ? Eigen::ColMajor : Eigen::RowMajor), max_rows, cols> &x)
    {
        assert(cols > 0);
        const int nrows = rows > 0 ? rows : x.size() / cols;
        Eigen::Matrix<double, rows, cols, (max_rows > 1 ? Eigen::ColMajor : Eigen::RowMajor), max_rows, cols> out;
        out.setZero(nrows, cols);
        
        for (int i = 0, id = 0; i < nrows; i++)
            for (int d = 0; d < cols; d++, id++)
                out(i, d) = x(id).getValue();

        return out;
    }

    enum class FD_RULE { CENTRAL, LEFT, RIGHT };
    
    void my_finite_gradient(const Eigen::VectorXd& x, const std::function<double(const Eigen::VectorXd&)> &f, Eigen::VectorXd &grad, FD_RULE rule = FD_RULE::CENTRAL, const double eps = 1e-7);

    // assume unit vector d
    double func1(
        const Eigen::Ref<const Eigen::Vector3d> &t,
        const Eigen::Ref<const Eigen::Vector3d> &d,
        const double &alpha, const double &beta);

    std::tuple<double, Vector6d> func1_grad(
        const Eigen::Ref<const Eigen::Vector3d> &t, 
        const Eigen::Ref<const Eigen::Vector3d> &d,
        const double &alpha, const double &beta);

    std::tuple<double, Vector6d, Matrix6d> func1_hess(
        const Eigen::Ref<const Eigen::Vector3d> &t, 
        const Eigen::Ref<const Eigen::Vector3d> &d,
        const double &alpha, const double &beta);
}
