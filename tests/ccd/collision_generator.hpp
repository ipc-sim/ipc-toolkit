#pragma once

#include <memory>

#include <Eigen/Core>
#include <catch2/catch.hpp>

struct TestImpact {
    Eigen::Vector2d p_t0, e0_t0, e1_t0;
    Eigen::Vector2d dp, de0, de1;
    double toi;
};

TestImpact generate_random_impact(const bool rigid);

class TestImpactGenerator : public Catch::Generators::IGenerator<TestImpact> {
protected:
    TestImpact current;
    bool rigid;
    size_t current_i, max_i;

public:
    TestImpactGenerator(size_t value, bool rigid);

    bool next() override;

    TestImpact const& get() const override;
};

Catch::Generators::GeneratorWrapper<TestImpact>
random_impacts(size_t value, bool rigid);
