#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_FILIB

#include <ipc/utils/eigen_ext.hpp>

#include <filib/interval.hpp>
#include <spdlog/fmt/bundled/ostream.h>

namespace filib {
/// @brief Interval type
class Interval : public interval {
public:
    /// @brief Construct an Interval from a filib interval
    /// @param i The filib interval
    Interval(const interval& i) : interval(i) { }

    /// @brief Construct an empty Interval
    Interval()
    {
        this->INF = 0;
        this->SUP = 0;
    }

    /// @brief Construct an Interval from a single value
    /// @param x The value
    explicit Interval(double x)
    {
        this->INF = x;
        this->SUP = x;
    }

    /// @brief Construct an Interval from two values
    /// @param x The infimum value
    /// @param y The supremum value
    Interval(double x, double y)
    {
        this->INF = x;
        this->SUP = y;
    }

    // friend std::ostream& operator<<(std::ostream& out, const Interval& i)
    // {
    //     return out << "[" << i.INF << ", " << i.SUP << "]";
    // }
};
} // namespace filib

template <> struct fmt::formatter<filib::Interval> : ostream_formatter { };

namespace ipc {

/// @brief 2D vector of intervals
using Vector2I = Eigen::Vector2<filib::Interval>;
/// @brief 3D vector of intervals
using Vector3I = Eigen::Vector3<filib::Interval>;
/// @brief Dynamic vector of intervals
using VectorXI = Eigen::VectorX<filib::Interval>;

/// @brief 2D row vector of intervals
using RowVector2I = Eigen::RowVector2<filib::Interval>;
/// @brief 3D row vector of intervals
using RowVector3I = Eigen::RowVector3<filib::Interval>;
/// @brief Dynamic row vector of intervals
using RowVectorXI = Eigen::RowVectorX<filib::Interval>;

/// @brief 2x2 matrix of intervals
using Matrix2I = Eigen::Matrix2<filib::Interval>;
/// @brief 3x3 matrix of intervals
using Matrix3I = Eigen::Matrix3<filib::Interval>;
/// @brief Dynamic matrix of intervals
using MatrixXI = Eigen::MatrixX<filib::Interval>;

/// @brief Dynamic vector of intervals with a maximum size of 3
using VectorMax3I = VectorMax3<filib::Interval>;
/// @brief Dynamic row vector of intervals with a maximum size of 3
using RowVectorMax3I = RowVectorMax3<filib::Interval>;
/// @brief Dynamic matrix of intervals with a maximum size of 3x3
using MatrixMax3I = MatrixMax3<filib::Interval>;

/// @brief Compute the squared L2 norm of an n-dimensional interval
/// @note This should be used instead of the .squaredNorm() method of Eigen because it avoids negative values in intermediate computations.
/// @param v The n-dimensional interval
/// @return The squared L2 norm of the interval
filib::Interval squared_norm(const Eigen::Ref<const VectorXI>& v); // L2 norm

/// @brief Compute the L2 norm of a n-dimensional interval
/// @note This should be used instead of the .norm() method of Eigen because it avoids negative values in intermediate computations.
/// @param v The n-dimensional interval
/// @return The L2 norm of the interval
filib::Interval norm(const Eigen::Ref<const VectorXI>& v); // L2 norm

} // namespace ipc

namespace Eigen {

template <typename BinOp>
struct ScalarBinaryOpTraits<filib::Interval, double, BinOp> {
    typedef filib::Interval ReturnType;
};

template <typename BinOp>
struct ScalarBinaryOpTraits<double, filib::Interval, BinOp> {
    typedef filib::Interval ReturnType;
};
} // namespace Eigen

#endif