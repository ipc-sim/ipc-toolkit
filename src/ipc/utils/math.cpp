#include "math.hpp"
#include <math.h>

const static double alpha = 4.; // control smoothness

namespace ipc {
    template <typename scalar>
    scalar smooth_heaviside(const scalar &x)
    {
        return 1. / (1 + exp(-alpha * x));
    }

    template double smooth_heaviside(const double &x);

    template <typename scalar>
    scalar heaviside(const scalar &x)
    {
        if (x > 0)
            return 1;
        else
            return 0;
    }

    template double heaviside(const double &x);

    template <typename scalar>
    scalar smooth_max(const scalar &x, const scalar &y)
    {
        scalar a = exp(alpha * x);
        scalar b = exp(alpha * y);
        return (x * a + y * b) / (a + b);
    }

    template double smooth_max(const double &, const double &);

    template <typename scalar>
    scalar smooth_min(const scalar &x, const scalar &y)
    {
        scalar a = exp(-alpha * x);
        scalar b = exp(-alpha * y);
        return (x * a + y * b) / (a + b);
    }

    template double smooth_min(const double &, const double &);

    template <typename scalar>
    scalar cubic_spline(const scalar &x)
    {
        if (x <= -2)
            return 0;
        if (x <= -1)
            return (x + 2)*(x + 2)*(x + 2) / 6;
        if (x <= 0)
            return 2. / 3 - x * x * (1 + x / 2);
        if (x <= 1)
            return 2. / 3 - x * x * (1 - x / 2);
        if (x < 2)
            return -(x - 2)*(x - 2)*(x - 2) / 6;
        return 0;
    }

    template double cubic_spline(const double &);

    template <typename scalar>
    scalar inv_barrier(const scalar &x, const double eps, const double r)
    {
        return cubic_spline(2 * x / eps) / abs(pow(x, r));
    }

    template double inv_barrier(const double &, const double, const double);
}