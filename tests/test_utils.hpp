#pragma once

#include <string>

#include <Eigen/Core>

#include <ipc/collision_constraint.hpp>

// Flatten the matrix rowwise
Eigen::VectorXd flatten(const Eigen::MatrixXd& X);

/// Unflatten rowwise
Eigen::MatrixXd unflatten(const Eigen::VectorXd& x, int dim);

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
