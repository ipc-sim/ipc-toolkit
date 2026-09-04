#include <catch2/catch_all.hpp>

#include <finitediff.hpp>

#include <ipc/dynamics/to_affine.hpp>

#include <iostream>
#include <memory>

using namespace ipc;
using namespace ipc::dynamics;

namespace {
std::shared_ptr<ToAffine>
make_to_affine(const bool rigid_cov, const int dim, const size_t num_bodies)
{
    if (rigid_cov) {
        return std::make_shared<RigidToAffine>(dim, num_bodies);
    }
    return std::make_shared<AffineToAffine>(dim, num_bodies);
}
} // namespace

TEST_CASE("ToAffine chain rule", "[to_affine][gradient][hessian]")
{
    const int dim = GENERATE(2, 3);
    const bool rigid_cov = GENERATE(true, false);
    const size_t N = 2;

    CAPTURE(dim, rigid_cov);

    auto to_affine = make_to_affine(rigid_cov, dim, N);
    const int rndof = to_affine->reduced_ndof();
    const int andof = to_affine->affine_ndof();

    // Random reduced DOFs; entries in [-1, 1] keep rotations inside the log
    // map.
    const Eigen::VectorXd x = Eigen::VectorXd::Random(rndof * N);

    // A random block-diagonal SPD affine-space quadratic φ(y) = ½yᵀCy + b·y.
    std::vector<Eigen::Triplet<double>> C_triplets;
    for (size_t i = 0; i < N; ++i) {
        const Eigen::MatrixXd Mi = Eigen::MatrixXd::Random(andof, andof);
        const Eigen::MatrixXd Ci =
            Mi.transpose() * Mi + Eigen::MatrixXd::Identity(andof, andof);
        for (int r = 0; r < andof; ++r) {
            for (int c = 0; c < andof; ++c) {
                C_triplets.emplace_back(i * andof + r, i * andof + c, Ci(r, c));
            }
        }
    }
    Eigen::SparseMatrix<double> C(andof * N, andof * N);
    C.setFromTriplets(C_triplets.begin(), C_triplets.end());
    const Eigen::VectorXd b = Eigen::VectorXd::Random(andof * N);

    const auto energy = [&](const Eigen::VectorXd& x_) -> double {
        const Eigen::VectorXd y = to_affine->to_affine(x_);
        return 0.5 * y.dot(C * y) + b.dot(y);
    };
    const auto grad_y = [&](const Eigen::VectorXd& y) -> Eigen::VectorXd {
        return C * y + b;
    };
    const auto analytic_grad =
        [&](const Eigen::VectorXd& x_) -> Eigen::VectorXd {
        return to_affine->apply_gradient(x_, grad_y(to_affine->to_affine(x_)));
    };

    // --- Round trip ---------------------------------------------------------
    if (rigid_cov) {
        CHECK(
            (to_affine->from_affine(to_affine->to_affine(x)) - x).norm()
            < 1e-9);
    } else {
        CHECK((to_affine->to_affine(x) - x).norm() == 0.0);
    }

    // --- Gradient: (dy/dx)ᵀ ∇φ vs FD of φ∘to_affine -------------------------
    const Eigen::VectorXd g = analytic_grad(x);
    Eigen::VectorXd fd_g;
    fd::finite_gradient(x, energy, fd_g);
    CHECK(fd::compare_gradient(g, fd_g));
    if (!fd::compare_gradient(g, fd_g)) {
        std::cout << "analytic:\n" << g.transpose() << "\n";
        std::cout << "numerical:\n" << fd_g.transpose() << "\n";
    }

    // --- Hessian: apply_hessian vs FD of the analytic gradient --------------
    const Eigen::SparseMatrix<double> H = to_affine->apply_hessian(
        x, grad_y(to_affine->to_affine(x)), C, PSDProjectionMethod::NONE);
    Eigen::MatrixXd fd_H;
    fd::finite_jacobian(x, analytic_grad, fd_H);
    CHECK(fd::compare_hessian(Eigen::MatrixXd(H), fd_H, 1e-3));
    if (!fd::compare_hessian(Eigen::MatrixXd(H), fd_H, 1e-3)) {
        std::cout << "analytic:\n" << Eigen::MatrixXd(H) << "\n\n";
        std::cout << "numerical:\n" << fd_H << "\n";
    }
}
