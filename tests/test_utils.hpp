#pragma once

#include <catch2/catch.hpp>

#include <string>

#include <Eigen/Core>

#include <ipc/collision_constraint.hpp>

#include <ipc/broad_phase/broad_phase.hpp>

#define GENERATE_BROAD_PHASE_METHODS()                                         \
    static_cast<BroadPhaseMethod>(                                             \
        GENERATE(range(0, static_cast<int>(BroadPhaseMethod::NUM_METHODS))));

static const std::string TEST_DATA_DIR(TEST_DATA_DIR_CSTR);

bool load_mesh(
    const std::string& mesh_name,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& E,
    Eigen::MatrixXi& F);

void mmcvids_to_constraints(
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXi& mmcvids,
    ipc::Constraints& constraints);

class RotationGenerator
    : public Catch::Generators::IGenerator<Eigen::Matrix3d> {
public:
    // Attempts to move the generator to the next element.
    // Returns true if successful (and thus has another element that can be
    // read)
    bool next() override;

    // Precondition:
    // The generator is either freshly constructed or the last call to next()
    // returned true
    Eigen::Matrix3d const& get() const override;

    static Catch::Generators::GeneratorWrapper<Eigen::Matrix3d> create();

protected:
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
};
