#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/potentials/friction_potential.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/utils/logger.hpp>

#include <igl/edges.h>

using namespace ipc;

void mmcvids_to_friction_collisions(
    Eigen::MatrixXi& E,
    Eigen::MatrixXi& F,
    const Eigen::MatrixXi& mmcvids,
    Eigen::VectorXd normal_force_magnitudes,
    const Eigen::MatrixXd& closest_points,
    const Eigen::MatrixXd& tangent_bases,
    TangentialCollisions& collisions)
{
    std::vector<Eigen::Vector2i> edges;
    std::vector<Eigen::Vector3i> faces;
    for (int i = 0; i < mmcvids.rows(); i++) {
        const auto mmcvid = mmcvids.row(i);
        TangentialCollision* collision;

        if (mmcvid[0] >= 0) { // Is EE?
            edges.emplace_back(mmcvid[0], mmcvid[1]);
            edges.emplace_back(mmcvid[2], mmcvid[3]);
            collisions.ee_collisions.emplace_back(
                edges.size() - 2, edges.size() - 1);
            collision = &(collisions.ee_collisions.back());
        } else {
            if (mmcvid[2] < 0) { // Is VV?
                collisions.vv_collisions.emplace_back(
                    -mmcvid[0] - 1, mmcvid[1]);
                CHECK(-mmcvid[3] >= 1);
                collisions.vv_collisions.back().weight = -mmcvid[3];
                normal_force_magnitudes[i] /= -mmcvid[3];
                collision = &(collisions.vv_collisions.back());

            } else if (mmcvid[3] < 0) { // Is EV?
                edges.emplace_back(mmcvid[1], mmcvid[2]);
                collisions.ev_collisions.emplace_back(
                    edges.size() - 1, -mmcvid[0] - 1);
                CHECK(-mmcvid[3] >= 1);
                collisions.ev_collisions.back().weight = -mmcvid[3];
                normal_force_magnitudes[i] /= -mmcvid[3];
                collision = &(collisions.ev_collisions.back());

            } else { // Is FV.
                faces.emplace_back(mmcvid[1], mmcvid[2], mmcvid[3]);
                collisions.fv_collisions.emplace_back(
                    faces.size() - 1, -mmcvid[0] - 1);
                collision = &(collisions.fv_collisions.back());
            }
        }

        collision->closest_point = closest_points.row(i);
        collision->tangent_basis = tangent_bases.middleRows(3 * i, 3);
        collision->normal_force_magnitude = normal_force_magnitudes[i];
    }

    E.resize(edges.size(), 2);
    for (int i = 0; i < edges.size(); i++) {
        E.row(i) = edges[i];
    }
    F.resize(faces.size(), 3);
    for (int i = 0; i < faces.size(); i++) {
        F.row(i) = faces[i];
    }
}

bool read_ipc_friction_data(
    const std::string& filename,
    Eigen::MatrixXd& V_start,
    Eigen::MatrixXd& V_lagged,
    Eigen::MatrixXd& V_end,
    Eigen::MatrixXi& E,
    Eigen::MatrixXi& F,
    NormalCollisions& collisions,
    TangentialCollisions& tangential_collisions,
    double& dhat,
    double& barrier_stiffness,
    double& epsv_times_h,
    double& mu,
    double& potential,
    Eigen::VectorXd& grad,
    Eigen::SparseMatrix<double>& hess)
{
    nlohmann::json data;

    std::ifstream input(filename);
    if (input.good()) {
        data = nlohmann::json::parse(input, nullptr, false);
    } else {
        logger().error("Unable to open IPC friction data file: {}", filename);
        return false;
    }

    if (data.is_discarded()) {
        logger().error("IPC friction data JSON is invalid: {}", filename);
        return false;
    }

    // Parameters
    dhat = sqrt(data["dhat_squared"].get<double>());
    barrier_stiffness = data["barrier_stiffness"];
    epsv_times_h = sqrt(data["epsv_times_h_squared"].get<double>());
    mu = data["mu"];

    // Dissipative potential value
    potential = data["energy"];

    // Potential gradient
    grad = data["gradient"];

    // Potential hessian
    std::vector<Eigen::Triplet<double>> hessian_triplets;
    hessian_triplets.reserve(data["hessian_triplets"].size());
    for (std::tuple<long, long, double> triplet : data["hessian_triplets"]) {
        hessian_triplets.emplace_back(
            std::get<0>(triplet), std::get<1>(triplet), std::get<2>(triplet));
    }
    hess.resize(grad.size(), grad.size());
    hess.setFromTriplets(hessian_triplets.begin(), hessian_triplets.end());

    // Mesh
    V_start = data["V_start"];
    V_lagged = data["V_lagged"];
    V_end = data["V_end"];

    // MMVCIDs
    const Eigen::MatrixXi mmcvids = data["mmcvids"];

    const Eigen::VectorXd lambda = data["normal_force_magnitudes"];

    const Eigen::MatrixXd coords = data["closest_point_coordinates"];

    const Eigen::MatrixXd bases = data["tangent_bases"];

    mmcvids_to_friction_collisions(
        E, F, mmcvids, lambda, coords, bases, tangential_collisions);
    tests::mmcvids_to_collisions(E, F, mmcvids, collisions);
    for (int i = 0; i < tangential_collisions.size(); i++) {
        tangential_collisions[i].mu = mu;
    }

    return true;
}

TEST_CASE(
    "Compare IPC friction derivatives", "[friction][gradient][hessian][data]")
{
    Eigen::MatrixXd V_start, V_lagged, V_end;
    Eigen::MatrixXi E, F;
    NormalCollisions collisions;
    TangentialCollisions expected_friction_collisions;
    double dhat, barrier_stiffness, epsv_times_h, mu;
    double expected_potential;
    Eigen::VectorXd expected_grad;
    Eigen::SparseMatrix<double> expected_hess;

    std::string scene_folder;
    int file_number = 0;
    SECTION("cube_cube")
    {
        scene_folder = "friction/cube_cube";
        file_number = GENERATE(range(0, 446));
    }
    // SECTION("chain")
    // {
    //     scene_folder = "friction/chain";
    //     file_number = GENERATE(range(0, 401));
    // }
    CAPTURE(scene_folder, file_number);

    bool success = read_ipc_friction_data(
        (tests::DATA_DIR / scene_folder
         / fmt::format("friction_data_{:d}.json", file_number))
            .string(),
        V_start, V_lagged, V_end, E, F, collisions,
        expected_friction_collisions, dhat, barrier_stiffness, epsv_times_h, mu,
        expected_potential, expected_grad, expected_hess);
    REQUIRE(success);

    Eigen::MatrixXi face_edges;
    igl::edges(F, face_edges);
    E.conservativeResize(E.rows() + face_edges.rows(), E.cols());
    E.bottomRows(face_edges.rows()) = face_edges;
    CollisionMesh mesh(V_start, E, F);

    TangentialCollisions tangential_collisions;
    tangential_collisions.build(
        mesh, V_lagged, collisions, BarrierPotential(dhat), barrier_stiffness,
        mu);
    REQUIRE(tangential_collisions.size() == collisions.size());
    REQUIRE(
        tangential_collisions.size() == expected_friction_collisions.size());

    const FrictionPotential D(epsv_times_h);

    REQUIRE(V_start.size() == V_lagged.size());
    REQUIRE(V_start.size() == V_end.size());
    REQUIRE(V_start.size() == expected_grad.size());

    for (int i = 0; i < tangential_collisions.size(); i++) {
        CAPTURE(i);
        const TangentialCollision& collision = tangential_collisions[i];
        const TangentialCollision& expected_collision =
            expected_friction_collisions[i];
        if (collision.closest_point.size() == 1) {
            CHECK(
                collision.closest_point[0]
                == Catch::Approx(expected_collision.closest_point[0]));
        } else {
            CHECK(collision.closest_point.isApprox(
                expected_collision.closest_point, 1e-12));
        }
        CHECK(collision.tangent_basis.isApprox(
            expected_collision.tangent_basis, 1e-12));
        CHECK(
            collision.normal_force_magnitude
            == Catch::Approx(expected_collision.normal_force_magnitude));
        CHECK(collision.mu == Catch::Approx(expected_collision.mu));

        CHECK(
            collision.vertex_ids(E, F) == expected_collision.vertex_ids(E, F));
    }

    const Eigen::MatrixXd velocity = V_end - V_start;

    double potential = D(tangential_collisions, mesh, velocity);

    CHECK(potential == Catch::Approx(expected_potential));

    Eigen::VectorXd grad = D.gradient(tangential_collisions, mesh, velocity);

    CHECK(grad.isApprox(expected_grad));

    Eigen::SparseMatrix<double> hess =
        D.hessian(tangential_collisions, mesh, velocity);

    CHECK(hess.isApprox(expected_hess));
}
