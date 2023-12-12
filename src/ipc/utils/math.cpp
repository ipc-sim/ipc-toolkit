#include "math.hpp"
#include <math.h>

const static double alpha = 4.; // control smoothness

namespace ipc {
    template <class scalar>
    scalar smooth_heaviside(const scalar &x)
    {
        return 1. / (1 + exp(-alpha * x));
    }

    template double smooth_heaviside<double>(const double &x);
    template float smooth_heaviside<float>(const float &x);

    template <class scalar>
    scalar heaviside(const scalar &x)
    {
        if (x > 0)
            return 1;
        else
            return 0;
    }

    template double heaviside<double>(const double &x);
    template float heaviside<float>(const float &x);

    template <class scalar>
    scalar smooth_max(const scalar &x, const scalar &y)
    {
        scalar a = exp(alpha * x);
        scalar b = exp(alpha * y);
        return (x * a + y * b) / (a + b);
    }

    template double smooth_max<double>(const double &, const double &);
    template float smooth_max<float>(const float &, const float &);

    template <class scalar>
    scalar smooth_min(const scalar &x, const scalar &y)
    {
        scalar a = exp(-alpha * x);
        scalar b = exp(-alpha * y);
        return (x * a + y * b) / (a + b);
    }

    template double smooth_min<double>(const double &, const double &);
    template float smooth_min<float>(const float &, const float &);
}