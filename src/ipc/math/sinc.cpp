#include "sinc.hpp"

namespace ipc {

namespace {
    // We use these bounds because for example 1 + x^2 = 1 for x < sqrt(ϵ).
    constexpr double TAYLOR_0_BOUND = std::numeric_limits<double>::epsilon();
    const double TAYLOR_2_BOUND = sqrt(TAYLOR_0_BOUND);
    const double TAYLOR_N_BOUND = sqrt(TAYLOR_2_BOUND);

    // WARNING: Assumes x is a single value and uses interval arithmetic to
    // account for rounding.
    filib::Interval sinc_interval_taylor(double x_double)
    {
        const filib::Interval x(x_double);

        if (abs(x.INF) >= TAYLOR_N_BOUND) {
            return sin(x) / x;
        }

        // approximation by taylor series in x at 0 up to order 5
        // 1 - x² / 6 + x⁴ / 120 = 1 + x² / 6 * (x² / 20 - 1)
        const filib::Interval squared_x = sqr(x);
        return 1.0 + squared_x / 6.0 * (squared_x / 20.0 - 1.0);
    }

    // Compute sinc'(x) / x
    inline double dsinc_over_x(double x)
    {
        static const double eps = 1e-4;

        double x2 = x * x;
        if (abs(x) > eps) {
            return (x * cos(x) - sin(x)) / (x2 * x);
        }

        // approximation by taylor series in x at 0 up to order 5
        return x2 * (-x2 / 840.0 + 1.0 / 30.0) - 1.0 / 3.0;
    }

    // Compute sinc"(x) / x² - sinc'(x) / x³
    inline double ddsinc_over_x2_minus_dsinc_over_x3(double x)
    {
        static const double eps = 0.1;

        double x2 = x * x;
        double x4 = x2 * x2;
        if (abs(x) > eps) {
            return ((3 - x2) * sin(x) - 3 * x * cos(x)) / (x4 * x);
        }

        // approximation by taylor series in x at 0 up to order 5
        return x4 / 7560.0 - x2 / 210.0 + 1.0 / 15.0;
    }
} // namespace

double sinc(const double x)
{
    if (abs(x) >= TAYLOR_N_BOUND) {
        return sin(x) / x;
    }

    // approximation by taylor series in x at 0 up to order 1
    double result = 1;

    if (abs(x) >= TAYLOR_0_BOUND) {
        const double squared_x = x * x;

        // approximation by taylor series in x at 0 up to order 3
        result -= squared_x / 6.0;

        if (abs(x) >= TAYLOR_2_BOUND) {
            // approximation by taylor series in x at 0 up to order 5
            result += (squared_x * squared_x) / 120.0;
        }
    }

    return result;
}

filib::Interval sinc(const filib::Interval& x)
{
    // Define two regions and use even symmetry of sinc.

    // A bound on sinc where it is monotonic ([0, ~4.4934])
    constexpr double monotonic_bound = 4.4934094579;

    // A conservative lower bound for sinc(x)
    // https://www.wolframalpha.com/input/?i=min+sin%28x%29%2Fx+between+4+and+5
    const filib::Interval bounds(-0.217233628211221659, 1.0);

    filib::Interval y = filib::Interval::empty(), x_pos = x;
    if (x.INF < 0) {
        if (x.SUP <= 0) {
            return sinc(-x); // sinc is an even function
        }
        // Split
        y = sinc(filib::Interval(0, -x.INF));
        x_pos = filib::Interval(0, x.SUP);
    }

    // Split the domain into two intervals:
    //  1) x ∩ [0, monotonic_bound]
    //  2) x ∩ [monotonic_bound, ∞)

    // Case 1 (Monotonic):
    filib::Interval x_gt_monotonic = x_pos;
    if (x_pos.INF <= monotonic_bound) {
        filib::Interval x_monotonic = x_pos;
        if (x_monotonic.SUP > monotonic_bound) {
            x_monotonic = filib::Interval(x_monotonic.INF, monotonic_bound);
            x_gt_monotonic = filib::Interval(monotonic_bound, x_pos.SUP);
        } else {
            x_gt_monotonic = filib::Interval::empty();
        }

        // sinc is monotonically decreasing, so flip lower() and upper().
        filib::Interval monotonic_y = (sinc_interval_taylor(x_monotonic.SUP)
                                       | sinc_interval_taylor(x_monotonic.INF))
            & bounds;

        // Combine the two intervals as a convex hull
        y |= monotonic_y;
    }

    // Case 2 (Not necessarily monotonic):
    if (!empty(x_gt_monotonic)) {
        // x_gt_monotonic is larger than one, so the division should be well
        // behaved.
        y |= sin(x_gt_monotonic) / x_gt_monotonic;
    }

    return intsec(y, bounds);
}

VectorMax3d sinc_norm_x_grad(Eigen::ConstRef<VectorMax3d> x)
{
    return dsinc_over_x(x.norm()) * x;
}

MatrixMax3d sinc_norm_x_hess(Eigen::ConstRef<VectorMax3d> x)
{
    double norm_x = x.norm();
    return ddsinc_over_x2_minus_dsinc_over_x3(norm_x) * x * x.transpose()
        + dsinc_over_x(norm_x) * MatrixMax3d::Identity(x.size(), x.size());
}

} // namespace ipc
