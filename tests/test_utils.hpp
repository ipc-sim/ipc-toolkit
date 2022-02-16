#pragma once

#include <catch2/catch.hpp>

#include <string>

#include <Eigen/Core>

#include <ipc/collision_constraint.hpp>

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
    virtual bool next();

    // Precondition:
    // The generator is either freshly constructed or the last call to next()
    // returned true
    virtual Eigen::Matrix3d const& get() const;

    static Catch::Generators::GeneratorWrapper<Eigen::Matrix3d> create();

protected:
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
};
