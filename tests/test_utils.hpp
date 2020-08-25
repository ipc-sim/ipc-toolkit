#pragma once

#include <string>

#include <Eigen/Core>

// Flatten the matrix rowwise
Eigen::VectorXd flatten(const Eigen::MatrixXd& X);

/// Unflatten rowwise
Eigen::MatrixXd unflatten(const Eigen::VectorXd& x, int dim);

bool load_mesh(
    const std::string& mesh_name,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& E,
    Eigen::MatrixXi& F);
