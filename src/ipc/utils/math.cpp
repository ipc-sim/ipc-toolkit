#include "math.hpp"
#include <math.h>
#include "AutodiffTypes.hpp"

// const static double alpha = 4.; // control smoothness

namespace ipc {
    // template <typename scalar>
    // scalar smooth_heaviside(const scalar &x)
    // {
    //     return 1. / (1 + exp(-alpha * x));
    // }

    // template double smooth_heaviside(const double &x);
    // // template AutodiffScalarGrad smooth_heaviside(const AutodiffScalarGrad &x);
    // // template AutodiffScalarHessian smooth_heaviside(const AutodiffScalarHessian &x);

    // template <typename scalar>
    // scalar heaviside(const scalar &x)
    // {
    //     if (x > 0)
    //         return scalar(1);
    //     else
    //         return scalar(0);
    // }

    // template double heaviside(const double &x);
    // // template AutodiffScalarGrad heaviside(const AutodiffScalarGrad &x);
    // // template AutodiffScalarHessian heaviside(const AutodiffScalarHessian &x);

    // template <typename scalar>
    // scalar smooth_max(const scalar &x, const scalar &y)
    // {
    //     scalar a = exp(alpha * x);
    //     scalar b = exp(alpha * y);
    //     return (x * a + y * b) / (a + b);
    // }

    // template double smooth_max(const double &, const double &);
    // // template AutodiffScalarGrad smooth_max(const AutodiffScalarGrad &, const AutodiffScalarGrad &);
    // // template AutodiffScalarHessian smooth_max(const AutodiffScalarHessian &, const AutodiffScalarHessian &);

    // template <typename scalar>
    // scalar smooth_min(const scalar &x, const scalar &y)
    // {
    //     scalar a = exp(-alpha * x);
    //     scalar b = exp(-alpha * y);
    //     return (x * a + y * b) / (a + b);
    // }

    // template double smooth_min(const double &, const double &);
    // // template AutodiffScalarGrad smooth_min(const AutodiffScalarGrad &, const AutodiffScalarGrad &);
    // // template AutodiffScalarHessian smooth_min(const AutodiffScalarHessian &, const AutodiffScalarHessian &);

    template <typename scalar>
    scalar cubic_spline(const scalar &x)
    {
        if (x <= -2)
            return scalar(0.);
        if (x <= -1)
            return (x + 2)*(x + 2)*(x + 2) / 6;
        if (x <= 0)
            return 2. / 3 - x * x * (1 + x / 2);
        if (x <= 1)
            return 2. / 3 - x * x * (1 - x / 2);
        if (x < 2)
            return -(x - 2)*(x - 2)*(x - 2) / 6;
        return scalar(0.);
    }

    template double cubic_spline(const double &);
    template AutodiffScalarGrad cubic_spline(const AutodiffScalarGrad &);
    template AutodiffScalarHessian cubic_spline(const AutodiffScalarHessian &);

    template <typename scalar>
    scalar inv_barrier(const scalar &x, const double eps, const double r)
    {
        return cubic_spline(2 * x / eps) / pow(sqrt(x*x), r);
    }

    template double inv_barrier(const double &, const double, const double);
    template AutodiffScalarGrad inv_barrier(const AutodiffScalarGrad &, const double, const double);
    template AutodiffScalarHessian inv_barrier(const AutodiffScalarHessian &, const double, const double);

}