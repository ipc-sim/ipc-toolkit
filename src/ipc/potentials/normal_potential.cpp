#include "normal_potential.hpp"

#include <ipc/utils/local_to_global.hpp>

#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

namespace ipc {

// -- Cumulative methods -------------------------------------------------------

Eigen::VectorXd NormalPotential::gauss_newton_hessian_diagonal(
    const NormalCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (collisions.empty()) {
        return Eigen::VectorXd::Zero(vertices.size());
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();
    const int dim = vertices.cols();

    tbb::combinable<Eigen::VectorXd> diag_storage(
        Eigen::VectorXd::Zero(vertices.size()));

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const NormalCollision& collision = collisions[i];

                const VectorMax12d local_diag =
                    this->gauss_newton_hessian_diagonal(
                        collision, collision.dof(vertices, edges, faces));

                const auto vids = collision.vertex_ids(edges, faces);

                // Don't be confused by the "gradient" in the name -- this just
                // scatters a local vector into a global vector.
                local_gradient_to_global_gradient(
                    local_diag, vids, dim, diag_storage.local());
            }
        });

    return diag_storage.combine([](const Eigen::VectorXd& a,
                                   const Eigen::VectorXd& b) { return a + b; });
}

double NormalPotential::gauss_newton_hessian_quadratic_form(
    const NormalCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    Eigen::ConstRef<Eigen::VectorXd> p) const
{
    assert(vertices.rows() == mesh.num_vertices());
    assert(p.size() == vertices.size());

    if (collisions.empty()) {
        return 0.0;
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();
    const int dim = vertices.cols();

    // Reshape p into a matrix for convenient DOF extraction
    const auto p_mat =
        p.reshaped<VERTEX_DERIVATIVE_LAYOUT>(vertices.rows(), dim);

    return tbb::parallel_reduce(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()), 0.0,
        [&](const tbb::blocked_range<size_t>& r, double partial_sum) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                partial_sum += this->gauss_newton_hessian_quadratic_form(
                    collisions[i], collisions[i].dof(vertices, edges, faces),
                    collisions[i].dof(p_mat, edges, faces));
            }
            return partial_sum;
        },
        std::plus<double>());
}

Eigen::SparseMatrix<double> NormalPotential::shape_derivative(
    const NormalCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (collisions.empty()) {
        return Eigen::SparseMatrix<double>(vertices.size(), vertices.size());
    }

    const Eigen::MatrixXd& rest_positions = mesh.rest_positions();
    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    const int ndof = vertices.size();

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                this->shape_derivative(
                    collisions[i], collisions[i].vertex_ids(edges, faces),
                    collisions[i].dof(rest_positions, edges, faces),
                    collisions[i].dof(vertices, edges, faces), local_triplets,
                    mesh.num_vertices());
            }
        });

    Eigen::SparseMatrix<double> shape_derivative(ndof, ndof);
    for (const auto& local_triplets : storage) {
        Eigen::SparseMatrix<double> local_shape_derivative(ndof, ndof);
        local_shape_derivative.setFromTriplets(
            local_triplets.begin(), local_triplets.end());
        shape_derivative += local_shape_derivative;
    }
    return shape_derivative;
}

// -- Single collision methods -------------------------------------------------

double NormalPotential::operator()(
    const NormalCollision& collision,
    Eigen::ConstRef<VectorMax12d> positions) const
{
    // w * m(x) * f(d(x))
    // NOTE: can save a multiplication by checking if !collision.is_mollified()
    const double d = collision.compute_distance(positions);
    return collision.weight * collision.mollifier(positions)
        * (*this)(d, collision.dmin);
}

VectorMax12d NormalPotential::gradient(
    const NormalCollision& collision,
    Eigen::ConstRef<VectorMax12d> positions) const
{
    // d(x)
    const double d = collision.compute_distance(positions);

    // Defensive guard: for mollified collisions (e.g. edge-edge near-parallel),
    // d=0 means the mollifier m=0 and ∇m=0, so the full gradient
    // f(d)·∇m + m·f'(d)·∇d = 0·∇m + 0·f'(0)·∇d = 0 regardless of f'(0).
    // Evaluating f'(0) directly would assert-fail inside BarrierPotential, so
    // we short-circuit here. d=0 should never occur after the
    // PARALLEL_THRESHOLD fix in edge_edge_distance_type, but this guards
    // against future regressions.
    if (collision.is_mollified() && d <= 0) {
        return VectorMax12d::Zero(positions.size());
    }

    // ∇d(x)
    const VectorMax12d grad_d = collision.compute_distance_gradient(positions);

    // f(d(x))
    const double f = (*this)(d, collision.dmin);
    // f'(d(x))
    const double grad_f = gradient(d, collision.dmin);

    if (!collision.is_mollified()) {
        // ∇[f(d(x))] = f'(d(x)) * ∇d(x)
        return (collision.weight * grad_f) * grad_d;
    }

    const double m = collision.mollifier(positions); // m(x)
    const VectorMax12d grad_m =
        collision.mollifier_gradient(positions); // ∇m(x)

    // ∇[m(x) * f(d(x))] = f(d(x)) * ∇m(x) + m(x) * ∇ f(d(x))
    return (collision.weight * f) * grad_m
        + (collision.weight * m * grad_f) * grad_d;
}

MatrixMax12d NormalPotential::hessian(
    const NormalCollision& collision,
    Eigen::ConstRef<VectorMax12d> positions,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    // d(x)
    const double d = collision.compute_distance(positions);

    // Defensive guard: same reasoning as in gradient() — for mollified
    // collisions at d=0, every term contains a factor of m=0 or ∇m=0, so the
    // full hessian is zero. Evaluating f'(0) or f''(0) would crash.
    if (collision.is_mollified() && d <= 0) {
        return MatrixMax12d::Zero(positions.size(), positions.size());
    }

    // ∇d(x)
    const VectorMax12d grad_d = collision.compute_distance_gradient(positions);
    // ∇²d(x)
    const MatrixMax12d hess_d = collision.compute_distance_hessian(positions);

    // f'(d(x))
    const double grad_f = gradient(d, collision.dmin);
    // f"(d(x))
    const double hess_f = hessian(d, collision.dmin);

    MatrixMax12d hess;
    if (!collision.is_mollified()) {
        // ∇²[f(d(x))] = ∇(f'(d(x)) * ∇d(x))
        //             = f"(d(x)) * ∇d(x) * ∇d(x)ᵀ + f'(d(x)) * ∇²d(x)
        hess = (collision.weight * hess_f) * grad_d * grad_d.transpose()
            + (collision.weight * grad_f) * hess_d;
    } else {
        const double f = (*this)(d, collision.dmin); // f(d(x))

        // m(x)
        const double m = collision.mollifier(positions);
        // ∇ m(x)
        const VectorMax12d grad_m = collision.mollifier_gradient(positions);
        // ∇² m(x)
        const MatrixMax12d hess_m = collision.mollifier_hessian(positions);

        const double weighted_m = collision.weight * m;

        // ∇f(d(x)) * ∇m(x)ᵀ
        const MatrixMax12d grad_f_grad_m =
            (collision.weight * grad_f) * grad_d * grad_m.transpose();

        // ∇²[m(x) * f(d(x))] = ∇[∇m(x) * f(d(x)) + m(x) * ∇f(d(x))]
        //                    = ∇²m(x) * f(d(x)) + ∇f(d(x)) * ∇m(x)ᵀ
        //                      + ∇m(x) * ∇f(d(x))ᵀ + m(x) * ∇²f(d(x))
        hess = (collision.weight * f) * hess_m + grad_f_grad_m
            + grad_f_grad_m.transpose()
            + (weighted_m * hess_f) * grad_d * grad_d.transpose()
            + (weighted_m * grad_f) * hess_d;
    }

    // Need to project entire hessian because w can be negative
    return project_to_psd(hess, project_hessian_to_psd);
}

// -- Single collision distance-vector methods ---------------------------------

VectorMax12d NormalPotential::gauss_newton_hessian_diagonal(
    const NormalCollision& collision,
    Eigen::ConstRef<VectorMax12d> positions) const
{
    const int d = collision.dim(positions.size());

    // Compute coefficients and distance vector together
    VectorMax4d coeffs;
    const VectorMax3d t = collision.compute_distance_vector(positions, coeffs);

    // d² = tᵀt
    const double d_sqr = t.squaredNorm();
    // f'(d²)
    const double grad_f = gradient(d_sqr, collision.dmin);
    // f''(d²)
    const double hess_f = hessian(d_sqr, collision.dmin);

    // diag(∂²b/∂x²) ≈ w · [4·f"(d²)·diag((∂t/∂x·t)(∂t/∂x·t)ᵀ)
    //                       + 2·f'(d²)·diag((∂t/∂x)(∂t/∂x)ᵀ)]
    //
    // From Eq. 11: diag((∂t/∂x)(∂t/∂x)ᵀ) = [c₀²,...,c₀²,c₁²,...,cₙ²]
    //   (each cᵢ² repeated dim times)
    // From Eq. 12: diag((∂t/∂x·t)(∂t/∂x·t)ᵀ) = [c₀t, c₁t, ..., cₙt]^∘2
    //   (element-wise square)

    const VectorMax12d diag_JtJ =
        CollisionStencil::diag_distance_vector_outer(coeffs, d);
    const VectorMax12d diag_JtttJt =
        CollisionStencil::diag_distance_vector_t_outer(coeffs, t);

    const double s1 = collision.weight * 4.0 * hess_f;
    const double s2 = collision.weight * 2.0 * grad_f;

    return s1 * diag_JtttJt + s2 * diag_JtJ;
}

double NormalPotential::gauss_newton_hessian_quadratic_form(
    const NormalCollision& collision,
    Eigen::ConstRef<VectorMax12d> positions,
    Eigen::ConstRef<VectorMax12d> p) const
{
    const int d = collision.dim(positions.size());

    // Compute coefficients and distance vector together
    VectorMax4d coeffs;
    const VectorMax3d t = collision.compute_distance_vector(positions, coeffs);

    // d² = tᵀt
    const double d_sqr = t.squaredNorm();
    // f'(d²)
    const double grad_f = gradient(d_sqr, collision.dmin);
    // f''(d²)
    const double hess_f = hessian(d_sqr, collision.dmin);

    // pᵀHp ≈ w · [4·f"(d²)·(pᵀ·(∂t/∂x)·t)² + 2·f'(d²)·‖pᵀ·(∂t/∂x)‖²]
    //
    // From Eq. 13: pᵀ(∂t/∂x·t)(∂t/∂x·t)ᵀp = (pᵀ·∂t/∂x·t)²
    // From Eq. 14: pᵀ(∂t/∂x)(∂t/∂x)ᵀp = ‖pᵀ·∂t/∂x‖²
    //
    // pᵀ·∂t/∂x = ∑ cᵢ pᵢ  (a dim-dimensional vector)

    const VectorMax3d pT_J =
        CollisionStencil::contract_distance_vector_jacobian(coeffs, p, d);

    // (pᵀ·∂t/∂x·t)² = (pT_J · t)²
    const double pTJt = pT_J.dot(t);
    const double term1 = pTJt * pTJt;

    // ‖pᵀ·∂t/∂x‖² = ‖pT_J‖²
    const double term2 = pT_J.squaredNorm();

    return collision.weight * (4.0 * hess_f * term1 + 2.0 * grad_f * term2);
}

void NormalPotential::shape_derivative(
    const NormalCollision& collision,
    const std::array<index_t, 4>& vertex_ids,
    Eigen::ConstRef<VectorMax12d> rest_positions, // = x̄
    Eigen::ConstRef<VectorMax12d> positions,      // = x̄ + u
    std::vector<Eigen::Triplet<double>>& out,
    const int n_total_verts) const
{
    assert(rest_positions.size() == positions.size());

    const int dim = collision.dim(positions.size());
    assert(positions.size() % collision.num_vertices() == 0);

    // Compute:
    // ∇ₓ (w ∇ᵤf(d(x̄+u))) = (∇ₓw)(∇ᵤf(d(x̄+u)))ᵀ + w ∇ₓ∇ᵤf(d(x̄+u))
    //                         (first term)        (second term)

    // First term:
    if (collision.weight_gradient.size() <= 0) {
        throw std::runtime_error(
            "Shape derivative is not computed for collisions!");
    }

    if (collision.weight_gradient.nonZeros() > 0) {
        VectorMax12d grad_f = gradient(collision, positions);
        assert(collision.weight != 0);
        grad_f.array() /= collision.weight; // remove weight

        for (int i = 0; i < collision.num_vertices(); i++) {
            for (int d = 0; d < dim; d++) {
                int row_idx;
                if constexpr (VERTEX_DERIVATIVE_LAYOUT == Eigen::RowMajor) {
                    row_idx = vertex_ids[i] * dim + d;
                } else {
                    assert(n_total_verts > 0);
                    row_idx = n_total_verts * d + vertex_ids[i];
                }
                using Itr = Eigen::SparseVector<double>::InnerIterator;
                for (Itr j(collision.weight_gradient); j; ++j) {
                    out.emplace_back(
                        row_idx, j.index(), grad_f[dim * i + d] * j.value());
                }
            }
        }
    }

    // Second term:
    MatrixMax12d local_hess;
    if (!collision.is_mollified()) {
        // w ∇ₓ∇ᵤf = w ∇ᵤ²f
        local_hess = hessian(collision, positions, PSDProjectionMethod::NONE);
    } else {
        // d(x̄+u)
        const double d = collision.compute_distance(positions);
        // ∇d(x̄+u)
        const VectorMax12d grad_d =
            collision.compute_distance_gradient(positions);
        // ∇²d(x̄+u)
        const MatrixMax12d hess_d =
            collision.compute_distance_hessian(positions);

        // f(d(x̄+u))
        const double f = collision.weight * (*this)(d, collision.dmin);
        // ∇ᵤ f(d(x̄+u))
        const Vector12d gradu_f =
            (collision.weight * gradient(d, collision.dmin)) * grad_d;
        // ∇ᵤ² f(d(x̄+u))
        const Matrix12d hessu_f =
            (collision.weight * hessian(d, collision.dmin)) * grad_d
                * grad_d.transpose()
            + (collision.weight * gradient(d, collision.dmin)) * hess_d;

        // ε(x̄)
        const double eps_x = collision.mollifier_threshold(rest_positions);

        // m(x̄,u)
        const double m = collision.mollifier(positions, eps_x);
        // ∇ᵤ m(x̄,u)
        const Vector12d gradu_m =
            collision.mollifier_gradient(positions, eps_x);
        // ∇ₓ m(x̄,u)
        const Vector12d gradx_m =
            collision.mollifier_gradient_wrt_x(rest_positions, positions);
        // ∇ₓ∇ᵤ m(x̄,u)
        const Matrix12d jac_m = collision.mollifier_gradient_jacobian_wrt_x(
            rest_positions, positions);

        // Only compute the second term of the shape derivative
        // ∇ₓ (f ∇ᵤm + m ∇ᵤf) = f ∇ₓ∇ᵤm + (∇ₓm)(∇ᵤf)ᵀ + (∇ᵤf)(∇ᵤm)ᵀ + m ∇ᵤ²f
        local_hess = f * jac_m + gradx_m * gradu_f.transpose()
            + gradu_f * gradu_m.transpose() + m * hessu_f;
    }

    local_hessian_to_global_triplets(
        local_hess, vertex_ids, dim, out, n_total_verts);
}

} // namespace ipc