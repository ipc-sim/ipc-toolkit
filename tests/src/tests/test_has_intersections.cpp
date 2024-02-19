#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>

#include <ipc/ipc.hpp>

#include <igl/edges.h>

using namespace ipc;

Eigen::MatrixXi remove_faces_with_degenerate_edges(
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
    Eigen::MatrixXi F_new(F.rows(), F.cols());
    int num_faces = 0;
    for (int i = 0; i < F.rows(); ++i) {
        if (V.row(F(i, 0)) != V.row(F(i, 1)) && V.row(F(i, 0)) != V.row(F(i, 2))
            && V.row(F(i, 1)) != V.row(F(i, 2))) {
            F_new.row(num_faces++) = F.row(i);
        }
    }
    F_new.conservativeResize(num_faces, F.cols());
    return F_new;
}

bool combine_meshes(
    const std::string& mesh1_name,
    const std::string& mesh2_name,
    const Eigen::Matrix3d& R1,
    const Eigen::Matrix3d& R2,
    int dim,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& E,
    Eigen::MatrixXi& F)
{
    Eigen::MatrixXd V1, V2;
    Eigen::MatrixXi E1, E2, F1, F2;

    bool success = tests::load_mesh(mesh1_name, V1, E1, F1)
        && tests::load_mesh(mesh2_name, V2, E2, F2);
    if (!success) {
        return false;
    }

    V1 = V1 * R1.transpose(); // (RVᵀ)ᵀ = VRᵀ
    V2 = V2 * R2.transpose();

    REQUIRE(dim <= V1.cols());
    REQUIRE(dim <= V2.cols());
    V = Eigen::MatrixXd(V1.rows() + V2.rows(), dim);
    V.topRows(V1.rows()) = V1.leftCols(dim);
    V.bottomRows(V2.rows()) = V2.leftCols(dim);

    REQUIRE(E1.cols() == E2.cols());
    E = Eigen::MatrixXi(E1.rows() + E2.rows(), E1.cols());
    E.topRows(E1.rows()) = E1;
    E.bottomRows(E2.rows()) = E2;
    E.bottomRows(E2.rows()).array() += V1.rows();

    REQUIRE(F1.cols() == F2.cols());
    F = Eigen::MatrixXi(F1.rows() + F2.rows(), F1.cols());
    F.topRows(F1.rows()) = F1;
    F.bottomRows(F2.rows()) = F2;
    F.bottomRows(F2.rows()).array() += V1.rows();

    F = remove_faces_with_degenerate_edges(V, F);
    if (F.rows() != 0) {
        igl::edges(F, E);
    }

    return true;
}

TEST_CASE("Has intersections", "[intersection]")
{
    std::string mesh1_name = GENERATE("cube.obj", "bunny.obj");
    std::string mesh2_name = GENERATE("cube.obj", "bunny.obj");
    int dim = GENERATE(2, 3);

#ifdef NDEBUG
    Eigen::Matrix3d R1 = GENERATE(take(4, tests::RotationGenerator::create()));
    Eigen::Matrix3d R2 = GENERATE(take(4, tests::RotationGenerator::create()));
#else
    Eigen::Matrix3d R1 = GENERATE(take(2, tests::RotationGenerator::create()));
    Eigen::Matrix3d R2 = GENERATE(take(2, tests::RotationGenerator::create()));
#endif

    const BroadPhaseMethod broad_phase_method = GENERATE_BROAD_PHASE_METHODS();

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    bool success = combine_meshes(mesh1_name, mesh2_name, R1, R2, dim, V, E, F);
    REQUIRE(success);

    CAPTURE(broad_phase_method);
    CHECK(has_intersections(CollisionMesh(V, E, F), V, broad_phase_method));
}
