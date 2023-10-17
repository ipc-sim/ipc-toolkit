#include <catch2/catch_test_macros.hpp>

#include <ipc/collision_mesh.hpp>

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

    CollisionMesh mesh(V, E, Eigen::MatrixXi(), W);

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
    Eigen::VectorXd gf = mesh.to_full_dof(g);
    Eigen::VectorXd expected_gf(6);
    expected_gf << 1, 1, -2, 0, 0, -2;
    CHECK(gf == expected_gf);

    Eigen::MatrixXd H(8, 8);
    H << 1, 0, 0, 0, 0, 0, 0, 0, //
        0, 1, 0, 0, 0, 0, 0, 0,  //
        0, 0, 1, 0, 0, 0, 0, 0,  //
        0, 0, 0, 1, 0, 0, 0, 0,  //
        0, 0, 0, 0, 1, 0, 0, 0,  //
        0, 0, 0, 0, 0, 1, 0, 0,  //
        0, 0, 0, 0, 0, 0, 1, 0,  //
        0, 0, 0, 0, 0, 0, 0, 1;  //
    Eigen::MatrixXd Hf =
        mesh.to_full_dof(Eigen::SparseMatrix<double>(H.sparseView()));

    Eigen::MatrixXd expected_Hf(6, 6);
    expected_Hf << 1, 0, 0, 0, 0, 0, //
        0, 1, 0, 0, 0, 0,            //
        0, 0, 2, 0, 1, 0,            //
        0, 0, 0, 2, 0, 1,            //
        0, 0, 1, 0, 2, 0,            //
        0, 0, 0, 1, 0, 2;            //

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

    CollisionMesh mesh(V, Eigen::MatrixXi(), Eigen::MatrixXi());

    Eigen::VectorXi expected_codim_vertices(4);
    expected_codim_vertices << 0, 1, 2, 3;
    CHECK(mesh.codim_vertices() == expected_codim_vertices);
}