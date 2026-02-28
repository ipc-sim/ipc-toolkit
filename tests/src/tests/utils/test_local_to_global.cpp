#include <ipc/utils/local_to_global.hpp>

#include <catch2/catch_test_macros.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>

using namespace ipc;

// ─── helpers ────────────────────────────────────────────────────────────────

/// Convert a RowMajor-ordered dense vector to ColMajor ordering.
///   RowMajor:  [x0 y0 z0  x1 y1 z1  …]
///   ColMajor:  [x0 x1 …  y0 y1 …  z0 z1 …]
static Eigen::VectorXd
row_to_col(const Eigen::VectorXd& v, int n_verts, int dim)
{
    return v.reshaped<Eigen::RowMajor>(n_verts, dim)
        .reshaped<Eigen::ColMajor>();
}

/// Convert a RowMajor-ordered sparse matrix to ColMajor DOF ordering.
///   RowMajor rows/cols indexed as (dim*v + d)
///   ColMajor rows/cols indexed as (n_verts*d + v)
static Eigen::MatrixXd
row_to_col_mat(const Eigen::MatrixXd& M, int n_verts, int dim)
{
    const int N = n_verts * dim;
    assert(M.rows() == N && M.cols() == N);
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(N, N);
    for (int vi = 0; vi < n_verts; ++vi) {
        for (int di = 0; di < dim; ++di) {
            for (int vj = 0; vj < n_verts; ++vj) {
                for (int dj = 0; dj < dim; ++dj) {
                    out(n_verts * di + vi, n_verts * dj + vj) =
                        M(dim * vi + di, dim * vj + dj);
                }
            }
        }
    }
    return out;
}

/// Same as row_to_col_mat but for a non-square Jacobian where row and column
/// vertex counts may differ.
static Eigen::MatrixXd row_to_col_jac(
    const Eigen::MatrixXd& M, int n_row_verts, int n_col_verts, int dim)
{
    const int Nr = n_row_verts * dim;
    const int Nc = n_col_verts * dim;
    assert(M.rows() == Nr && M.cols() == Nc);
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(Nr, Nc);
    for (int vi = 0; vi < n_row_verts; ++vi) {
        for (int di = 0; di < dim; ++di) {
            for (int vj = 0; vj < n_col_verts; ++vj) {
                for (int dj = 0; dj < dim; ++dj) {
                    out(n_row_verts * di + vi, n_col_verts * dj + vj) =
                        M(dim * vi + di, dim * vj + dj);
                }
            }
        }
    }
    return out;
}

/// Assemble a dense matrix from a vector of triplets.
static Eigen::MatrixXd triplets_to_dense(
    const std::vector<Eigen::Triplet<double>>& trips, int rows, int cols)
{
    Eigen::SparseMatrix<double> sp(rows, cols);
    sp.setFromTriplets(trips.begin(), trips.end());
    return Eigen::MatrixXd(sp);
}

// ─── Test parameters ────────────────────────────────────────────────────────

// Global mesh: 5 vertices in 3D  → 15 DOFs.
// Local element touches vertices {1, 3} (2 verts, 6 local DOFs).
static constexpr int DIM = 3;
static constexpr int N_TOTAL_VERTS = 5;
static constexpr int N_DOF = N_TOTAL_VERTS * DIM;

static const std::vector<long> ELEM_IDS = { { 1, 3 } };
static constexpr int N_LOCAL_VERTS = 2;
static constexpr int N_LOCAL_DOF = N_LOCAL_VERTS * DIM;

/// Deterministic local gradient (6 entries).
static Eigen::VectorXd make_local_grad()
{
    Eigen::VectorXd g(N_LOCAL_DOF);
    // RowMajor local ordering: [x1 y1 z1 x3 y3 z3]
    g << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
    return g;
}

/// Deterministic local hessian (6×6).
static Eigen::MatrixXd make_local_hessian()
{
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(N_LOCAL_DOF, N_LOCAL_DOF);
    for (int i = 0; i < N_LOCAL_DOF; ++i) {
        for (int j = 0; j < N_LOCAL_DOF; ++j) {
            H(i, j) = 10.0 * (i + 1) + (j + 1); // nonzero everywhere
        }
    }
    return H;
}

// ═══════════════════════════════════════════════════════════════════════════
//  1. local_gradient_to_global_gradient – dense Eigen::VectorXd
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE(
    "local_gradient_to_global_gradient dense RowMajor vs ColMajor",
    "[utils][local_to_global][gradient][dense]")
{
    const Eigen::VectorXd local_grad = make_local_grad();

    // --- RowMajor assembly ---
    Eigen::VectorXd grad_row = Eigen::VectorXd::Zero(N_DOF);
    local_gradient_to_global_gradient<std::vector<long>, Eigen::RowMajor>(
        local_grad, ELEM_IDS, DIM, grad_row);

    // Verify a few values manually.
    // Vertex 1, dim 0 → index 3*1+0 = 3
    CHECK(grad_row[3] == 1.0);
    // Vertex 1, dim 1 → index 3*1+1 = 4
    CHECK(grad_row[4] == 2.0);
    // Vertex 3, dim 2 → index 3*3+2 = 11
    CHECK(grad_row[11] == 6.0);

    // --- ColMajor assembly ---
    Eigen::VectorXd grad_col = Eigen::VectorXd::Zero(N_DOF);
    local_gradient_to_global_gradient<std::vector<long>, Eigen::ColMajor>(
        local_grad, ELEM_IDS, DIM, grad_col);

    // ColMajor: vertex v, dim d → index N_TOTAL_VERTS*d + v
    // Vertex 1, dim 0 → index 5*0+1 = 1
    CHECK(grad_col[1] == 1.0);
    // Vertex 1, dim 1 → index 5*1+1 = 6
    CHECK(grad_col[6] == 2.0);
    // Vertex 3, dim 2 → index 5*2+3 = 13
    CHECK(grad_col[13] == 6.0);

    // Structural check: both vectors should be the permutation of each other.
    Eigen::VectorXd expected_col = row_to_col(grad_row, N_TOTAL_VERTS, DIM);
    CHECK((grad_col.array() == expected_col.array()).all());
}

TEST_CASE(
    "local_gradient_to_global_gradient dense accumulates from two elements",
    "[utils][local_to_global][gradient][dense]")
{
    // Two elements that share vertex 3.
    const std::vector<long> ids_a = { 0, 3 };
    const std::vector<long> ids_b = { 3, 4 };

    Eigen::VectorXd local_a(N_LOCAL_DOF), local_b(N_LOCAL_DOF);
    local_a << 1, 2, 3, 10, 20, 30;
    local_b << 100, 200, 300, 7, 8, 9;

    Eigen::VectorXd grad_row = Eigen::VectorXd::Zero(N_DOF);
    local_gradient_to_global_gradient<std::vector<long>, Eigen::RowMajor>(
        local_a, ids_a, DIM, grad_row);
    local_gradient_to_global_gradient<std::vector<long>, Eigen::RowMajor>(
        local_b, ids_b, DIM, grad_row);

    Eigen::VectorXd grad_col = Eigen::VectorXd::Zero(N_DOF);
    local_gradient_to_global_gradient<std::vector<long>, Eigen::ColMajor>(
        local_a, ids_a, DIM, grad_col);
    local_gradient_to_global_gradient<std::vector<long>, Eigen::ColMajor>(
        local_b, ids_b, DIM, grad_col);

    // Vertex 3 receives contributions from both elements.
    // RowMajor vertex 3 dim 0 → index 9: 10 + 100 = 110
    CHECK(grad_row[9] == 110.0);
    // ColMajor vertex 3 dim 0 → index 5*0+3 = 3: should also be 110
    CHECK(grad_col[3] == 110.0);

    Eigen::VectorXd expected_col = row_to_col(grad_row, N_TOTAL_VERTS, DIM);
    CHECK((grad_col.array() == expected_col.array()).all());
}

// ═══════════════════════════════════════════════════════════════════════════
//  2. local_gradient_to_global_gradient – sparse Eigen::SparseVector
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE(
    "local_gradient_to_global_gradient sparse RowMajor vs ColMajor",
    "[utils][local_to_global][gradient][sparse]")
{
    const Eigen::VectorXd local_grad = make_local_grad();

    Eigen::SparseVector<double> grad_row(N_DOF);
    local_gradient_to_global_gradient<std::vector<long>, Eigen::RowMajor>(
        local_grad, ELEM_IDS, DIM, grad_row);

    Eigen::SparseVector<double> grad_col(N_DOF);
    local_gradient_to_global_gradient<std::vector<long>, Eigen::ColMajor>(
        local_grad, ELEM_IDS, DIM, grad_col);

    // Convert both to dense for easy comparison.
    Eigen::VectorXd dr = Eigen::VectorXd(grad_row);
    Eigen::VectorXd dc = Eigen::VectorXd(grad_col);

    Eigen::VectorXd expected_col = row_to_col(dr, N_TOTAL_VERTS, DIM);
    CHECK((dc.array() == expected_col.array()).all());

    // Spot-check: only 2*DIM entries should be nonzero.
    CHECK(grad_row.nonZeros() == N_LOCAL_DOF);
    CHECK(grad_col.nonZeros() == N_LOCAL_DOF);
}

// ═══════════════════════════════════════════════════════════════════════════
//  3. local_hessian_to_global_triplets – std::vector<Triplet>
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE(
    "local_hessian_to_global_triplets (triplet vector) RowMajor vs ColMajor",
    "[utils][local_to_global][hessian][triplet]")
{
    const Eigen::MatrixXd local_hess = make_local_hessian();

    // --- RowMajor ---
    std::vector<Eigen::Triplet<double>> trips_row;
    local_hessian_to_global_triplets<std::vector<long>, Eigen::RowMajor>(
        local_hess, ELEM_IDS, DIM, trips_row);

    // --- ColMajor ---
    std::vector<Eigen::Triplet<double>> trips_col;
    local_hessian_to_global_triplets<std::vector<long>, Eigen::ColMajor>(
        local_hess, ELEM_IDS, DIM, trips_col, N_TOTAL_VERTS);

    // Both should produce N_LOCAL_DOF^2 triplets.
    CHECK(trips_row.size() == N_LOCAL_DOF * N_LOCAL_DOF);
    CHECK(trips_col.size() == N_LOCAL_DOF * N_LOCAL_DOF);

    Eigen::MatrixXd H_row = triplets_to_dense(trips_row, N_DOF, N_DOF);
    Eigen::MatrixXd H_col = triplets_to_dense(trips_col, N_DOF, N_DOF);

    Eigen::MatrixXd expected_col = row_to_col_mat(H_row, N_TOTAL_VERTS, DIM);
    CHECK((H_col.array() == expected_col.array()).all());
}

TEST_CASE(
    "local_hessian_to_global_triplets (triplet vector) ColMajor spot-check",
    "[utils][local_to_global][hessian][triplet]")
{
    // 2 verts, dim=2, ids={0,2}, n_total=4
    const int dim = 2;
    const int n_total = 4;
    const std::vector<long> ids = { 0, 2 };

    // 4×4 local hessian
    Eigen::MatrixXd H(4, 4);
    H << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

    std::vector<Eigen::Triplet<double>> trips;
    local_hessian_to_global_triplets<std::vector<long>, Eigen::ColMajor>(
        H, ids, dim, trips, n_total);

    Eigen::MatrixXd M = triplets_to_dense(trips, n_total * dim, n_total * dim);

    // RowMajor local (i=0,d=0) → global ColMajor row = n_total*0 + ids[0] = 0
    // RowMajor local (i=0,d=1) → global ColMajor row = n_total*1 + ids[0] = 4
    // RowMajor local (i=1,d=0) → global ColMajor row = n_total*0 + ids[1] = 2
    // RowMajor local (i=1,d=1) → global ColMajor row = n_total*1 + ids[1] = 6

    // H(0,0) = local(dim*0+0, dim*0+0) = H(0,0)=1, row=0 col=0
    CHECK(M(0, 0) == 1.0);
    // H(0,1) = local(0,1) = 2, row=0, col=n_total*1+0=4
    CHECK(M(0, 4) == 2.0);
    // H(1,0) = local(1,0) = 5, row=n_total*1+0=4, col=0
    CHECK(M(4, 0) == 5.0);
    // H(2,3) = local(dim*1+0, dim*1+1) = H(2,3) = 12
    //   row = n_total*0+ids[1] = 2, col = n_total*1+ids[1] = 6
    CHECK(M(2, 6) == 12.0);
    // H(3,2) = local(dim*1+1, dim*1+0) = H(3,2) = 15
    //   row = n_total*1+ids[1] = 6, col = n_total*0+ids[1] = 2
    CHECK(M(6, 2) == 15.0);
}

// ═══════════════════════════════════════════════════════════════════════════
//  4. local_hessian_to_global_triplets – MatrixCache
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE(
    "local_hessian_to_global_triplets (MatrixCache) RowMajor vs ColMajor",
    "[utils][local_to_global][hessian][cache]")
{
    const Eigen::MatrixXd local_hess = make_local_hessian();

    // --- RowMajor via MatrixCache ---
    SparseMatrixCache cache_row(N_DOF, N_DOF);
    cache_row.reserve(N_LOCAL_DOF * N_LOCAL_DOF);
    local_hessian_to_global_triplets<
        Eigen::MatrixXd, std::vector<long>, Eigen::RowMajor>(
        local_hess, ELEM_IDS, DIM, cache_row);

    Eigen::MatrixXd H_row = Eigen::MatrixXd(cache_row.get_matrix());

    // --- ColMajor via MatrixCache ---
    SparseMatrixCache cache_col(N_DOF, N_DOF);
    cache_col.reserve(N_LOCAL_DOF * N_LOCAL_DOF);
    local_hessian_to_global_triplets<
        Eigen::MatrixXd, std::vector<long>, Eigen::ColMajor>(
        local_hess, ELEM_IDS, DIM, cache_col, N_TOTAL_VERTS);

    Eigen::MatrixXd H_col = Eigen::MatrixXd(cache_col.get_matrix());

    Eigen::MatrixXd expected_col = row_to_col_mat(H_row, N_TOTAL_VERTS, DIM);
    CHECK((H_col.array() == expected_col.array()).all());
}

TEST_CASE(
    "local_hessian_to_global_triplets (MatrixCache) ColMajor with two elements",
    "[utils][local_to_global][hessian][cache]")
{
    // Two elements that share vertex 2.
    const int dim = 2;
    const int n_total = 4;
    const int n_dof = n_total * dim;
    const std::vector<long> ids_a = { 0, 2 };
    const std::vector<long> ids_b = { 2, 3 };

    Eigen::MatrixXd Ha = Eigen::MatrixXd::Ones(4, 4) * 2.0;
    Eigen::MatrixXd Hb = Eigen::MatrixXd::Ones(4, 4) * 3.0;

    // RowMajor
    SparseMatrixCache cr(n_dof, n_dof);
    cr.reserve(32);
    local_hessian_to_global_triplets<
        Eigen::MatrixXd, std::vector<long>, Eigen::RowMajor>(
        Ha, ids_a, dim, cr);
    local_hessian_to_global_triplets<
        Eigen::MatrixXd, std::vector<long>, Eigen::RowMajor>(
        Hb, ids_b, dim, cr);
    Eigen::MatrixXd M_row = Eigen::MatrixXd(cr.get_matrix());

    // ColMajor
    SparseMatrixCache cc(n_dof, n_dof);
    cc.reserve(32);
    local_hessian_to_global_triplets<
        Eigen::MatrixXd, std::vector<long>, Eigen::ColMajor>(
        Ha, ids_a, dim, cc, n_total);
    local_hessian_to_global_triplets<
        Eigen::MatrixXd, std::vector<long>, Eigen::ColMajor>(
        Hb, ids_b, dim, cc, n_total);
    Eigen::MatrixXd M_col = Eigen::MatrixXd(cc.get_matrix());

    Eigen::MatrixXd expected = row_to_col_mat(M_row, n_total, dim);
    CHECK((M_col.array() == expected.array()).all());

    // The (2,2) block in RowMajor should accumulate 2+3=5 for the
    // diagonal dim entries of vertex 2.
    CHECK(M_row(dim * 2, dim * 2) == 5.0);
    CHECK(M_col(n_total * 0 + 2, n_total * 0 + 2) == 5.0);
}

// ═══════════════════════════════════════════════════════════════════════════
//  5. local_jacobian_to_global_triplets
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE(
    "local_jacobian_to_global_triplets RowMajor vs ColMajor square",
    "[utils][local_to_global][jacobian]")
{
    // Square case (same as hessian but using the jacobian API).
    const Eigen::MatrixXd local_jac = make_local_hessian();

    std::vector<Eigen::Triplet<double>> trips_row;
    local_jacobian_to_global_triplets<
        Eigen::MatrixXd, std::vector<long>, std::vector<long>, Eigen::RowMajor>(
        local_jac, ELEM_IDS, ELEM_IDS, DIM, trips_row);

    std::vector<Eigen::Triplet<double>> trips_col;
    local_jacobian_to_global_triplets<
        Eigen::MatrixXd, std::vector<long>, std::vector<long>, Eigen::ColMajor>(
        local_jac, ELEM_IDS, ELEM_IDS, DIM, trips_col, N_TOTAL_VERTS,
        N_TOTAL_VERTS);

    Eigen::MatrixXd J_row = triplets_to_dense(trips_row, N_DOF, N_DOF);
    Eigen::MatrixXd J_col = triplets_to_dense(trips_col, N_DOF, N_DOF);

    Eigen::MatrixXd expected = row_to_col_mat(J_row, N_TOTAL_VERTS, DIM);
    CHECK((J_col.array() == expected.array()).all());
}

TEST_CASE(
    "local_jacobian_to_global_triplets ColMajor non-square",
    "[utils][local_to_global][jacobian]")
{
    // Non-square: 2 row-verts (ids {0,3}), 3 col-verts (ids {1,2,4}), dim=2.
    const int dim = 2;
    const int n_row_total = 5; // global row vertices
    const int n_col_total = 5; // global col vertices
    const std::vector<long> row_ids = { 0, 3 };
    const std::vector<long> col_ids = { 1, 2, 4 };
    const int n_row_local = 2;
    const int n_col_local = 3;
    const int lr = n_row_local * dim; // 4
    const int lc = n_col_local * dim; // 6

    Eigen::MatrixXd J(lr, lc);
    for (int i = 0; i < lr; ++i)
        for (int j = 0; j < lc; ++j)
            J(i, j) = 100.0 * (i + 1) + (j + 1);

    // RowMajor
    std::vector<Eigen::Triplet<double>> trips_row;
    local_jacobian_to_global_triplets<
        Eigen::MatrixXd, std::vector<long>, std::vector<long>, Eigen::RowMajor>(
        J, row_ids, col_ids, dim, trips_row);

    // ColMajor
    std::vector<Eigen::Triplet<double>> trips_col;
    local_jacobian_to_global_triplets<
        Eigen::MatrixXd, std::vector<long>, std::vector<long>, Eigen::ColMajor>(
        J, row_ids, col_ids, dim, trips_col, n_row_total, n_col_total);

    const int Nr = n_row_total * dim;
    const int Nc = n_col_total * dim;
    Eigen::MatrixXd J_row = triplets_to_dense(trips_row, Nr, Nc);
    Eigen::MatrixXd J_col = triplets_to_dense(trips_col, Nr, Nc);

    Eigen::MatrixXd expected =
        row_to_col_jac(J_row, n_row_total, n_col_total, dim);
    CHECK((J_col.array() == expected.array()).all());

    // Spot-check a single entry.
    // J(0,0) = 101.  RowMajor: row=dim*row_ids[0]+0=0, col=dim*col_ids[0]+0=2
    CHECK(J_row(0, 2) == 101.0);
    // ColMajor: row=n_row_total*0+row_ids[0]=0, col=n_col_total*0+col_ids[0]=1
    CHECK(J_col(0, 1) == 101.0);
}

// ═══════════════════════════════════════════════════════════════════════════
//  6. Edge cases
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("ColMajor gradient 2D", "[utils][local_to_global][gradient][2D]")
{
    // 2D problem: 4 vertices, dim=2.
    const int dim = 2;
    const int n_total = 4;
    const int n_dof = n_total * dim; // 8
    const std::vector<long> ids = { 1, 3 };
    Eigen::VectorXd local_grad(4);
    local_grad << 10, 20, 30, 40; // [x1 y1 x3 y3]

    // RowMajor
    Eigen::VectorXd gr = Eigen::VectorXd::Zero(n_dof);
    local_gradient_to_global_gradient<std::vector<long>, Eigen::RowMajor>(
        local_grad, ids, dim, gr);
    // Expected RowMajor: index 2*1+0=2 →10, 2*1+1=3 →20, 2*3+0=6 →30, 2*3+1=7
    // →40
    CHECK(gr[2] == 10.0);
    CHECK(gr[3] == 20.0);
    CHECK(gr[6] == 30.0);
    CHECK(gr[7] == 40.0);

    // ColMajor
    Eigen::VectorXd gc = Eigen::VectorXd::Zero(n_dof);
    local_gradient_to_global_gradient<std::vector<long>, Eigen::ColMajor>(
        local_grad, ids, dim, gc);
    // Expected ColMajor: n_total*0+1=1 →10, n_total*1+1=5 →20,
    //                    n_total*0+3=3 →30, n_total*1+3=7 →40
    CHECK(gc[1] == 10.0);
    CHECK(gc[5] == 20.0);
    CHECK(gc[3] == 30.0);
    CHECK(gc[7] == 40.0);

    Eigen::VectorXd expected = row_to_col(gr, n_total, dim);
    CHECK((gc.array() == expected.array()).all());
}

TEST_CASE(
    "ColMajor single vertex element",
    "[utils][local_to_global][gradient][single]")
{
    // A single vertex local element.
    const int dim = 3;
    const int n_total = 3;
    const int n_dof = n_total * dim;
    const std::vector<long> ids = { 2 };
    Eigen::VectorXd local_grad(3);
    local_grad << 7, 8, 9;

    Eigen::VectorXd gr = Eigen::VectorXd::Zero(n_dof);
    local_gradient_to_global_gradient<std::vector<long>, Eigen::RowMajor>(
        local_grad, ids, dim, gr);

    Eigen::VectorXd gc = Eigen::VectorXd::Zero(n_dof);
    local_gradient_to_global_gradient<std::vector<long>, Eigen::ColMajor>(
        local_grad, ids, dim, gc);

    Eigen::VectorXd expected = row_to_col(gr, n_total, dim);
    CHECK((gc.array() == expected.array()).all());

    // RowMajor: vertex 2, dim 0 → 3*2+0=6
    CHECK(gr[6] == 7.0);
    // ColMajor: vertex 2, dim 0 → 3*0+2=2
    CHECK(gc[2] == 7.0);
}

TEST_CASE(
    "ColMajor hessian identity-like local matrix",
    "[utils][local_to_global][hessian][identity]")
{
    // Local identity hessian to verify diagonal placement.
    const int dim = 2;
    const int n_total = 4;
    const int n_dof = n_total * dim;
    const std::vector<long> ids = { 0, 3 };

    Eigen::MatrixXd H = Eigen::MatrixXd::Identity(4, 4);

    std::vector<Eigen::Triplet<double>> trips_row;
    local_hessian_to_global_triplets<std::vector<long>, Eigen::RowMajor>(
        H, ids, dim, trips_row);

    std::vector<Eigen::Triplet<double>> trips_col;
    local_hessian_to_global_triplets<std::vector<long>, Eigen::ColMajor>(
        H, ids, dim, trips_col, n_total);

    Eigen::MatrixXd M_row = triplets_to_dense(trips_row, n_dof, n_dof);
    Eigen::MatrixXd M_col = triplets_to_dense(trips_col, n_dof, n_dof);

    // RowMajor diagonal entries at (0,0),(1,1),(6,6),(7,7).
    CHECK(M_row(0, 0) == 1.0);
    CHECK(M_row(1, 1) == 1.0);
    CHECK(M_row(6, 6) == 1.0);
    CHECK(M_row(7, 7) == 1.0);

    // ColMajor diagonal entries at (0,0),(4,4),(3,3),(7,7).
    CHECK(M_col(0, 0) == 1.0);
    CHECK(M_col(4, 4) == 1.0);
    CHECK(M_col(3, 3) == 1.0);
    CHECK(M_col(7, 7) == 1.0);

    Eigen::MatrixXd expected = row_to_col_mat(M_row, n_total, dim);
    CHECK((M_col.array() == expected.array()).all());
}

TEST_CASE(
    "ColMajor with extra ids in container",
    "[utils][local_to_global][gradient][extra_ids]")
{
    // ids container is larger than needed – only first n_verts are used.
    const int dim = 2;
    const int n_total = 5;
    const int n_dof = n_total * dim;
    const std::vector<long> ids = { 1, 4, 99 }; // 3 ids, but local only needs 2

    Eigen::VectorXd local_grad(4); // 2 verts * 2 dim
    local_grad << 5, 6, 7, 8;

    Eigen::VectorXd gr = Eigen::VectorXd::Zero(n_dof);
    local_gradient_to_global_gradient<std::vector<long>, Eigen::RowMajor>(
        local_grad, ids, dim, gr);

    Eigen::VectorXd gc = Eigen::VectorXd::Zero(n_dof);
    local_gradient_to_global_gradient<std::vector<long>, Eigen::ColMajor>(
        local_grad, ids, dim, gc);

    Eigen::VectorXd expected = row_to_col(gr, n_total, dim);
    CHECK((gc.array() == expected.array()).all());
}

TEST_CASE(
    "ColMajor with Eigen::VectorXi ids",
    "[utils][local_to_global][gradient][eigen_ids]")
{
    // Ensure IDContainer works with Eigen types.
    const int dim = 3;
    const int n_total = 4;
    const int n_dof = n_total * dim;
    Eigen::VectorXi ids(2);
    ids << 0, 3;

    Eigen::VectorXd local_grad(6);
    local_grad << 1, 2, 3, 4, 5, 6;

    Eigen::VectorXd gr = Eigen::VectorXd::Zero(n_dof);
    local_gradient_to_global_gradient<Eigen::VectorXi, Eigen::RowMajor>(
        local_grad, ids, dim, gr);

    Eigen::VectorXd gc = Eigen::VectorXd::Zero(n_dof);
    local_gradient_to_global_gradient<Eigen::VectorXi, Eigen::ColMajor>(
        local_grad, ids, dim, gc);

    Eigen::VectorXd expected = row_to_col(gr, n_total, dim);
    CHECK((gc.array() == expected.array()).all());
}
