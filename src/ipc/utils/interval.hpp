#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_FILIB

#include <ipc/utils/eigen_ext.hpp>

#include <filib/interval.hpp>

#include <spdlog/fmt/bundled/ostream.h>

namespace filib {
class Interval : public interval {
public:
    Interval(const interval& i) : interval(i) { }

    Interval()
    {
        this->INF = 0;
        this->SUP = 0;
    }

    Interval(double x)
    {
        this->INF = x;
        this->SUP = x;
    }

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

typedef Vector2<filib::Interval> Vector2I;
typedef Vector3<filib::Interval> Vector3I;
typedef VectorMax3<filib::Interval> VectorMax3I;
typedef VectorX<filib::Interval> VectorXI;
typedef RowVector2<filib::Interval> RowVector2I;
typedef RowVector3<filib::Interval> RowVector3I;
typedef RowVectorMax3<filib::Interval> RowVectorMax3I;
typedef RowVectorX<filib::Interval> RowVectorXI;
typedef Matrix2<filib::Interval> Matrix2I;
typedef Matrix3<filib::Interval> Matrix3I;
typedef MatrixMax3<filib::Interval> MatrixMax3I;
typedef MatrixX<filib::Interval> MatrixXI;

/// @brief Compute the L2 norm of a 3-dimensional interval
/// @param v The 3-dimensional interval
/// @return The L2 norm of the interval
filib::Interval squared_norm(const Eigen::Ref<const VectorXI>& v); // L2 norm

/// @brief Compute the L2 norm of a 3-dimensional interval
/// @param v The 3-dimensional interval
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