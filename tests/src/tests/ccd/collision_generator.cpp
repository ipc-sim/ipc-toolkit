#include "collision_generator.hpp"

#include <cstdlib> /* rand */

TestImpactGenerator::TestImpactGenerator(size_t value, bool _rigid)
    : rigid(_rigid)
    , current_i(0)
    , max_i(value)
{
    current = generate_random_impact(rigid);
}

bool TestImpactGenerator::next()
{
    current = generate_random_impact(rigid);
    return ++current_i < max_i;
}

TestImpact const& TestImpactGenerator::get() const { return current; }

Catch::Generators::GeneratorWrapper<TestImpact>
random_impacts(size_t value, bool rigid)
{
    return Catch::Generators::GeneratorWrapper<TestImpact>(
        Catch::Detail::make_unique<TestImpactGenerator>(value, rigid));
}

TestImpact generate_random_impact(const bool rigid)
{
    static const double kEPSILON = 1E-8;

    const auto rand_vec2 = []() {
        static const double kSCALE = 10.0;
        return kSCALE * Eigen::Vector2d::Random();
    };

    TestImpact impact;

    // set first edge random positions keeping in larger than kEPSILON
    impact.e0_t0 = rand_vec2();
    impact.e1_t0 = rand_vec2();

    Eigen::Vector2d e_dir = impact.e1_t0 - impact.e0_t0;
    double e_len = e_dir.norm();
    if (e_len < kEPSILON) {
        impact.e1_t0 += (e_dir) / e_len * kEPSILON;
        assert((impact.e1_t0 - impact.e0_t0).norm() >= kEPSILON);
    }

    // set first edge random velocities keeping in larger than kEPSILON.
    // note that this doest NOT ensure its length is > 0 over the full
    // time.
    impact.de0 = rand_vec2();
    if (rigid) {
        impact.de1 = impact.de0;
    } else {
        impact.de1 = rand_vec2();
        e_dir = (impact.e1_t0 + impact.de1) - (impact.e0_t0 + impact.de0);
        e_len = e_dir.norm();
        if (e_len < kEPSILON) {
            impact.de1 += (e_dir) / e_len * kEPSILON;
            assert(
                ((impact.e1_t0 + impact.de1) - (impact.e0_t0 + impact.de0))
                    .norm()
                >= kEPSILON);
        }
    }

    // Generate a time and position for the impact.
    // How do we ensure it is the first time of impact?
    impact.toi = double(rand()) / RAND_MAX;
    Eigen::Vector2d e0_t0_toi = impact.e0_t0 + impact.de0 * impact.toi;
    Eigen::Vector2d e1_t0_toi = impact.e1_t0 + impact.de1 * impact.toi;
    double alpha = double(rand()) / RAND_MAX;
    Eigen::Vector2d pc = e0_t0_toi + alpha * (e1_t0_toi - e0_t0_toi);

    // pc = p_t0 + dp * toi → (pc - p_t0)/toi = dp
    impact.p_t0 = rand_vec2();
    impact.dp = (pc - impact.p_t0) / impact.toi;

    return impact;
}
