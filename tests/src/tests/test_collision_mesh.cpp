#include <tests/dof_layout.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/config.hpp>

using namespace ipc;

TEST_CASE("Collision mesh", "[collision_mesh]")
{
    Eigen::MatrixXd V(4, 2);
    V << 0, 0, 1, 0, 0, 1, 1, 1;
    Eigen::MatrixXi E(4, 2);
    E << 0, 1, 1, 3, 3, 2, 2, 0;

    std::vector<Eigen::Triplet<double>> weights;
    weights.emplace_back(0, 0, 1);
    weights.emplace_back(1, 1, 1);
    weights.emplace_back(2, 2, 1);
    weights.emplace_back(3, 1, 1);
    weights.emplace_back(3, 2, 1);

    Eigen::SparseMatrix<double> W(4, 3);
    W.setFromTriplets(weights.begin(), weights.end());

    CollisionMesh mesh(V, E, /*F=*/Eigen::MatrixXi(), W);

    CHECK(mesh.max_edge_length() == Catch::Approx(1.).margin(1e-15));
    CHECK(mesh.edge_length(0) == Catch::Approx(1.).margin(1e-15));

    Eigen::MatrixXd U(3, 2);
    U << 0, 0, 1, 1, 0, 0;

    Eigen::MatrixXd Uc = mesh.map_displacements(U);
    Eigen::MatrixXd expected_Uc(4, 2);
    expected_Uc << 0, 0, 1, 1, 0, 0, 1, 1;
    CHECK(Uc == expected_Uc);

    Eigen::MatrixXd Vc = mesh.displace_vertices(U);
    Eigen::MatrixXd expected_Vc(4, 2);
    expected_Vc << 0, 0, 2, 1, 0, 1, 2, 2;
    CHECK(Vc == expected_Vc);

    Eigen::VectorXd g(8);
    g << 1, 1, -1, 1, 1, -1, -1, -1;
    g = tests::reorder_gradient(g, mesh.num_vertices(), V.cols());
    Eigen::VectorXd gf = mesh.to_full_dof(g);
    Eigen::VectorXd expected_gf(6);
    expected_gf << 1, 1, -2, 0, 0, -2;
    expected_gf = tests::reorder_gradient(expected_gf, U.rows(), U.cols());
    CHECK(gf == expected_gf);

    Eigen::MatrixXd H = Eigen::MatrixXd::Identity(8, 8);
    Eigen::MatrixXd Hf =
        mesh.to_full_dof(Eigen::SparseMatrix<double>(H.sparseView()));

    Eigen::MatrixXd expected_Hf(6, 6);
    expected_Hf << 1, 0, 0, 0, 0, 0, //
        0, 1, 0, 0, 0, 0,            //
        0, 0, 2, 0, 1, 0,            //
        0, 0, 0, 2, 0, 1,            //
        0, 0, 1, 0, 2, 0,            //
        0, 0, 0, 1, 0, 2;            //
    expected_Hf = tests::reorder_hessian(expected_Hf, U.rows(), U.cols());

    CHECK(Hf == expected_Hf);
}

TEST_CASE("Faces to edges", "[collision_mesh][faces_to_edges]")
{
    Eigen::MatrixXi F(1, 3), E(3, 2), expected_F2E(1, 3);
    F << 0, 1, 2;

    SECTION("Works")
    {
        SECTION("In order")
        {
            E << 0, 1, 1, 2, 2, 0;
            expected_F2E << 0, 1, 2;
        }

        SECTION("Reverse order")
        {
            E << 2, 0, 2, 1, 1, 0;
            expected_F2E << 2, 1, 0;
        }

        SECTION("Shuffled")
        {
            E << 0, 1, 2, 0, 2, 1;
            expected_F2E << 0, 2, 1;
        }

        CHECK(CollisionMesh::construct_faces_to_edges(F, E) == expected_F2E);
    }
    SECTION("Shouldnt work")
    {
        E << 0, 1, 1, 2, 0, 3;
        try {
            CollisionMesh::construct_faces_to_edges(F, E);
            FAIL("Should have thrown");
        } catch (const std::runtime_error& e) {
            SUCCEED("Should have thrown");
            CHECK(e.what() == std::string("Unable to find edge!"));
        } catch (...) {
            FAIL("Uknown exception thrown");
        }
    }
}

TEST_CASE("Codim points collision mesh", "[collision_mesh]")
{
    Eigen::MatrixXd V(4, 2);
    V << 0, 0, 1, 0, 0, 1, 1, 1;

    CollisionMesh mesh(V);

    Eigen::VectorXi expected_codim_vertices(4);
    expected_codim_vertices << 0, 1, 2, 3;
    CHECK(mesh.codim_vertices() == expected_codim_vertices);
}

TEST_CASE(
    "vertex_matrix_to_dof_matrix",
    "[collision_mesh][vertex_matrix_to_dof_matrix]")
{
    // Build a small 2×3 vertex-level matrix:
    //
    //   M_V = [ 1  0  2 ]
    //         [ 0  3  0 ]
    //
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.emplace_back(0, 0, 1.0);
    triplets.emplace_back(0, 2, 2.0);
    triplets.emplace_back(1, 1, 3.0);

    Eigen::SparseMatrix<double> M_V(2, 3);
    M_V.setFromTriplets(triplets.begin(), triplets.end());

    const int n_rows = static_cast<int>(M_V.rows()); // 2
    const int n_cols = static_cast<int>(M_V.cols()); // 3

    const int dim = GENERATE(2, 3);
    CAPTURE(dim);

    Eigen::SparseMatrix<double> M_dof =
        CollisionMesh::vertex_matrix_to_dof_matrix(M_V, dim);

    CHECK(M_dof.rows() == n_rows * dim);
    CHECK(M_dof.cols() == n_cols * dim);

    // Check the sparsity pattern: for each non-zero (r, c, v) in M_V,
    // it should appear at (global_row(r, d), global_col(c, d)) for
    // every d in [0, dim).
    Eigen::MatrixXd M_dense = Eigen::MatrixXd(M_dof);
    Eigen::MatrixXd expected =
        Eigen::MatrixXd::Zero(n_rows * dim, n_cols * dim);

    using InnerIterator = Eigen::SparseMatrix<double>::InnerIterator;
    for (int k = 0; k < M_V.outerSize(); ++k) {
        for (InnerIterator it(M_V, k); it; ++it) {
            for (int d = 0; d < dim; d++) {
                int gr, gc;
                if constexpr (VERTEX_DERIVATIVE_LAYOUT == Eigen::RowMajor) {
                    gr = dim * it.row() + d;
                    gc = dim * it.col() + d;
                } else {
                    gr = n_rows * d + it.row();
                    gc = n_cols * d + it.col();
                }
                expected(gr, gc) = it.value();
            }
        }
    }

    CHECK(M_dense == expected);

    // Verify semantics: M_dof * flatten(V) == flatten(M_V * V)
    // V is (n_cols × dim) because M_V has n_cols columns, one per vertex.
    Eigen::MatrixXd V(n_cols, dim); // NOLINT(bugprone-argument-comment)
    V.setRandom();

    Eigen::VectorXd x_dof = V.reshaped<VERTEX_DERIVATIVE_LAYOUT>();
    Eigen::VectorXd expected_dof =
        (M_V * V).reshaped<VERTEX_DERIVATIVE_LAYOUT>();

    Eigen::VectorXd actual_dof = M_dof * x_dof;

    CHECK(actual_dof.isApprox(expected_dof));
}
TEST_CASE("face_normals", "[collision_mesh][face_normals]")
{
    // Two triangles of a unit square in the z = 0 plane, wound oppositely so
    // their normals point in opposite directions.
    Eigen::MatrixXd V(4, 3);
    V << 0, 0, 0, //
        1, 0, 0,  //
        1, 1, 0,  //
        0, 1, 0;  //

    Eigen::MatrixXi F(2, 3);
    F << 0, 1, 2, // CCW  ⇒ +z
        0, 2, 3;  // CCW  ⇒ +z

    Eigen::MatrixXi E(5, 2);
    E << 0, 1, 1, 2, 2, 0, 2, 3, 3, 0;

    const CollisionMesh mesh(V, E, F);

    SECTION("rest positions")
    {
        const std::vector<Eigen::Vector3d> N = mesh.face_normals(V);

        REQUIRE(N.size() == mesh.num_faces());
        for (const Eigen::Vector3d& n : N) {
            CHECK(n.isApprox(Eigen::Vector3d::UnitZ()));
            CHECK(n.norm() == Catch::Approx(1.0).margin(1e-15));
        }
    }

    SECTION("normals follow deformed positions")
    {
        // Rotate the mesh 90° about the x-axis: +z normal ⇒ −y normal.
        Eigen::MatrixXd V_rot = V;
        for (int i = 0; i < V.rows(); ++i) {
            V_rot.row(i) << V(i, 0), -V(i, 2), V(i, 1);
        }

        const std::vector<Eigen::Vector3d> N = mesh.face_normals(V_rot);

        REQUIRE(N.size() == mesh.num_faces());
        for (const Eigen::Vector3d& n : N) {
            CHECK(n.isApprox(-Eigen::Vector3d::UnitY()));
        }
    }

    SECTION("winding flips the normal")
    {
        Eigen::MatrixXi F_flipped(2, 3);
        F_flipped << 0, 2, 1, // reversed ⇒ −z
            0, 3, 2;          // reversed ⇒ −z

        const CollisionMesh flipped(V, E, F_flipped);

        const std::vector<Eigen::Vector3d> N = flipped.face_normals(V);

        REQUIRE(N.size() == flipped.num_faces());
        for (const Eigen::Vector3d& n : N) {
            CHECK(n.isApprox(-Eigen::Vector3d::UnitZ()));
        }
    }
}
