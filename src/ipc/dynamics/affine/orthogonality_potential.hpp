#pragma once

#include <ipc/dynamics/affine/affine_body.hpp>

namespace ipc::affine {

class OrthogonalityPotential {
public:
    OrthogonalityPotential(const double stiffness) : stiffness(stiffness) { }
    virtual ~OrthogonalityPotential() = default;

    // -- Cumulative methods ---------------------------------------------------

    double operator()(const std::vector<AffineBody>& bodies) const;

    Eigen::VectorXd gradient(const std::vector<AffineBody>& bodies) const;

    Eigen::SparseMatrix<double> hessian(
        const std::vector<AffineBody>& bodies,
        const bool project_hessian_to_psd = false) const;

    // -- Single body methods ---------------------------------------------

    double operator()(const AffineBody& body) const;

    VectorMax12d gradient(const AffineBody& body) const;

    MatrixMax12d hessian(const AffineBody& body) const;

    double stiffness;
};

} // namespace ipc::affine