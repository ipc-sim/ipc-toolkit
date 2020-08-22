#include "edge_vertex_ccd_2D.hpp"

// namespace ipc {
namespace ccd {

/// \brief Value used to decide when a coeff (a,b,c) is zero.
static const double TOI_COEFF_EPSILON = 1E-12;
/// \brief Value used when solving for alpha given a toi, to avoid division
/// by zero.
static const double ALPHA_DIVISION_EPSILON = 1E-8;

namespace autogen {

    inline void edge_vertex_2D_time_of_impact_coeff(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk,
        double& a,
        double& b,
        double& c)
    {
        a = -Ui[0] * Uj[1] + Ui[0] * Uk[1] + Ui[1] * Uj[0] - Ui[1] * Uk[0]
            - Uj[0] * Uk[1] + Uj[1] * Uk[0];
        b = -Ui[0] * Vj[1] + Ui[0] * Vk[1] + Ui[1] * Vj[0] - Ui[1] * Vk[0]
            + Uj[0] * Vi[1] - Uj[0] * Vk[1] - Uj[1] * Vi[0] + Uj[1] * Vk[0]
            - Uk[0] * Vi[1] + Uk[0] * Vj[1] + Uk[1] * Vi[0] - Uk[1] * Vj[0];
        c = -Vi[0] * Vj[1] + Vi[0] * Vk[1] + Vi[1] * Vj[0] - Vi[1] * Vk[0]
            - Vj[0] * Vk[1] + Vj[1] * Vk[0];

        double max_coeff = std::max(a, std::max(b, c));
        a /= max_coeff;
        b /= max_coeff;
        c /= max_coeff;
    }

} // namespace autogen

bool compute_edge_vertex_time_of_impact(
    const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj,
    const Eigen::Vector2d& Uk,
    double& toi,
    double& alpha)
{
    double a, b, c;
    autogen::edge_vertex_2D_time_of_impact_coeff(
        Vi, Vj, Vk, Ui, Uj, Uk, a, b, c);

    // At most we will have two solutions for the quadratic equation
    // we initialize them with an invalid value (i.e default failure)
    // our approach uses two different formulas that solve the quadratic
    // equation
    //     (1)         (-b - sqrt(b^2 - 4 * a * c) / (2 * a)
    //     (2)         -2*c / (|b|+ sqrt(b^2 - 4 * a * c))

    auto is_not_zero = [](const double x) -> bool {
        return abs(x) > TOI_COEFF_EPSILON;
    };
    double x1(-1), x2(-1);
    double radicand = b * b - 4 * a * c;
    bool a_not_zero = is_not_zero(a);
    bool b_not_zero = is_not_zero(b);
    bool c_not_zero = is_not_zero(c);

    if (radicand > 0) {
        double sqrt_rad = sqrt(radicand);
        if (b > 0) {
            x1 = -2 * c / (b + sqrt_rad);
            // if (a_not_zero) {
            x2 = (-b - sqrt_rad) / (2 * a);
            //}
        } else { // note x1 and x2 switched
            x2 = 2 * c / (-b + sqrt_rad);
            // if (a_not_zero) {
            x1 = (-b + sqrt_rad) / (2 * a);
            //}
        }
    } else if (radicand == 0 && a_not_zero) {
        x1 = (-b) / (2 * a);
    } else if (!a_not_zero && !b_not_zero && !c_not_zero) {
        // check for case a=b=c=0

        // we will use the following approximation:
        //      distance to the closest edge vertex divided by rel
        //      velocity
        Eigen::Vector2d n = (Vj - Vi).normalized().cast<double>();
        double dist_ik = (Vi - Vk).norm();
        double dist_jk = (Vj - Vk).norm();

        if (dist_ik < dist_jk) { // closest i
            x1 = double(dist_ik) / (Uk - Ui).dot(n);
        } else {
            x1 = double(dist_jk) / (Uk - Uj).dot(-n);
        }
    }
    auto check_solution = [Vi, Vj, Vk, Ui, Uj, Uk](const double& t, double& s) {
        return t >= 0 && t <= 1
            && temporal_parameterization_to_spatial(
                   Vi, Vj, Vk, Ui, Uj, Uk, t, s)
            && s >= 0 && s <= 1;
    };

    // now check which of the solutions are valid
    double alpha1(-1), alpha2(-1);
    bool x1_valid = check_solution(x1, alpha1);
    bool x2_valid = check_solution(x2, alpha2);
    if (x1_valid && x2_valid) {
        toi = std::min(x1, x2);
        alpha = x1 < x2 ? alpha1 : alpha2;
    } else if (x1_valid) {
        toi = x1;
        alpha = alpha1;
    } else {
        toi = x2;
        alpha = alpha2;
    }

    return x1_valid || x2_valid;
}

bool temporal_parameterization_to_spatial(
    const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj,
    const Eigen::Vector2d& Uk,
    const double toi,
    double& alpha)
{
    // solves for alpha
    // Vk + Uk * toi =
    //      (Vi + Ui * toi) + alpha ((Vj + Uj * toi) - (Vi + Ui * toi)).

    auto numerator_0 = Vk[0] - Vi[0] + (Uk[0] - Ui[0]) * toi;
    auto numerator_1 = Vk[1] - Vi[1] + (Uk[1] - Ui[1]) * toi;
    auto denominator_0 = (Vj[0] - Vi[0] + (Uj[0] - Ui[0]) * toi);
    auto denominator_1 = (Vj[1] - Vi[1] + (Uj[1] - Ui[1]) * toi);

    // we need to divide to obtain alpha but we only need one dimension
    if (abs(denominator_0) > ALPHA_DIVISION_EPSILON) {
        alpha = numerator_0 / denominator_0;
        return true;
    }
    if (abs(denominator_1) > ALPHA_DIVISION_EPSILON) {
        alpha = numerator_1 / denominator_1;
        return true;
    }
    return false;
}

} // namespace ccd
// } // namespace ipc
