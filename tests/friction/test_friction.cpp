#include <catch2/catch_all.hpp>

#include <array>
#include <vector>

#include <finitediff.hpp>
#include <igl/edges.h>

#include <ipc/ipc.hpp>
#include <ipc/friction/friction_constraints.hpp>
#include <ipc/utils/logger.hpp>
#include <ipc/config.hpp>

#include "friction_data_generator.hpp"
#include "../utils.hpp"

#include <unsupported/Eigen/SparseExtra>

using namespace ipc;

TEST_CASE("Friction gradient and hessian", "[friction][gradient][hessian]")
{
    FrictionData data = friction_data_generator();
    const auto& [V0, V1, E, F, collision_constraints, mu, epsv_times_h, dhat, barrier_stiffness] =
        data;

    const Eigen::MatrixXd U = V1 - V0;

    const CollisionMesh mesh(V0, E, F);

    FrictionConstraints friction_constraints;
    friction_constraints.build(
        mesh, V0, collision_constraints, dhat, barrier_stiffness, mu);

    const Eigen::VectorXd grad =
        friction_constraints.compute_potential_gradient(mesh, U, epsv_times_h);

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        const Eigen::MatrixXd fd_U = fd::unflatten(x, data.V1.cols()) - data.V0;
        return friction_constraints.compute_potential(
            mesh, fd_U, data.epsv_times_h);
    };
    Eigen::VectorXd fgrad;
    fd::finite_gradient(fd::flatten(V1), f, fgrad);
    CHECK(fd::compare_gradient(grad, fgrad));

    const Eigen::MatrixXd hess =
        friction_constraints.compute_potential_hessian(mesh, U, epsv_times_h);
    Eigen::MatrixXd fhess;
    fd::finite_hessian(fd::flatten(V1), f, fhess);
    CHECK(fd::compare_hessian(hess, fhess, 1e-3));
}

///////////////////////////////////////////////////////////////////////////////

void mmcvids_to_friction_constraints(
    Eigen::MatrixXi& E,
    Eigen::MatrixXi& F,
    const Eigen::MatrixXi& mmcvids,
    Eigen::VectorXd normal_force_magnitudes,
    const Eigen::MatrixXd& closest_points,
    const Eigen::MatrixXd& tangent_bases,
    FrictionConstraints& constraints)
{
    std::vector<Eigen::Vector2i> edges;
    std::vector<Eigen::Vector3i> faces;
    for (int i = 0; i < mmcvids.rows(); i++) {
        const auto mmcvid = mmcvids.row(i);
        FrictionConstraint* constraint;

        if (mmcvid[0] >= 0) { // Is EE?
            edges.emplace_back(mmcvid[0], mmcvid[1]);
            edges.emplace_back(mmcvid[2], mmcvid[3]);
            constraints.ee_constraints.emplace_back(
                edges.size() - 2, edges.size() - 1);
            constraint = &(constraints.ee_constraints.back());
        } else {
            if (mmcvid[2] < 0) { // Is VV?
                constraints.vv_constraints.emplace_back(
                    -mmcvid[0] - 1, mmcvid[1]);
                CHECK(-mmcvid[3] >= 1);
                constraints.vv_constraints.back().weight = -mmcvid[3];
                normal_force_magnitudes[i] /= -mmcvid[3];
                constraint = &(constraints.vv_constraints.back());

            } else if (mmcvid[3] < 0) { // Is EV?
                edges.emplace_back(mmcvid[1], mmcvid[2]);
                constraints.ev_constraints.emplace_back(
                    edges.size() - 1, -mmcvid[0] - 1);
                CHECK(-mmcvid[3] >= 1);
                constraints.ev_constraints.back().weight = -mmcvid[3];
                normal_force_magnitudes[i] /= -mmcvid[3];
                constraint = &(constraints.ev_constraints.back());

            } else { // Is FV.
                faces.emplace_back(mmcvid[1], mmcvid[2], mmcvid[3]);
                constraints.fv_constraints.emplace_back(
                    faces.size() - 1, -mmcvid[0] - 1);
                constraint = &(constraints.fv_constraints.back());
            }
        }

        constraint->closest_point = closest_points.row(i);
        constraint->tangent_basis = tangent_bases.middleRows(3 * i, 3);
        constraint->normal_force_magnitude = normal_force_magnitudes[i];
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
    CollisionConstraints& collision_constraints,
    FrictionConstraints& friction_constraints,
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
    from_json(data["gradient"], grad);

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
    from_json(data["V_start"], V_start);
    from_json(data["V_lagged"], V_lagged);
    from_json(data["V_end"], V_end);

    // MMVCIDs
    Eigen::MatrixXi mmcvids;
    from_json(data["mmcvids"], mmcvids);

    Eigen::VectorXd lambda;
    from_json(data["normal_force_magnitudes"], lambda);

    Eigen::MatrixXd coords;
    from_json(data["closest_point_coordinates"], coords);

    Eigen::MatrixXd bases;
    from_json(data["tangent_bases"], bases);

    mmcvids_to_friction_constraints(
        E, F, mmcvids, lambda, coords, bases, friction_constraints);
    mmcvids_to_constraints(E, F, mmcvids, collision_constraints);
    for (int i = 0; i < friction_constraints.size(); i++) {
        friction_constraints[i].mu = mu;
    }

    return true;
}

TEST_CASE(
    "Compare IPC friction derivatives", "[friction][gradient][hessian][data]")
{
    Eigen::MatrixXd V_start, V_lagged, V_end;
    Eigen::MatrixXi E, F;
    CollisionConstraints collision_constraints;
    FrictionConstraints expected_friction_constraints;
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
        fmt::format(
            "{}{}/friction_data_{:d}.json", TEST_DATA_DIR, scene_folder,
            file_number),
        V_start, V_lagged, V_end, E, F, collision_constraints,
        expected_friction_constraints, dhat, barrier_stiffness, epsv_times_h,
        mu, expected_potential, expected_grad, expected_hess);
    REQUIRE(success);

    Eigen::MatrixXi face_edges;
    igl::edges(F, face_edges);
    E.conservativeResize(E.rows() + face_edges.rows(), E.cols());
    E.bottomRows(face_edges.rows()) = face_edges;
    CollisionMesh mesh(V_start, E, F);

    FrictionConstraints friction_constraints;
    friction_constraints.build(
        mesh, V_lagged, collision_constraints, dhat, barrier_stiffness, mu);

    REQUIRE(friction_constraints.size() == collision_constraints.size());
    REQUIRE(
        friction_constraints.size() == expected_friction_constraints.size());

    REQUIRE(V_start.size() == V_lagged.size());
    REQUIRE(V_start.size() == V_end.size());
    REQUIRE(V_start.size() == expected_grad.size());

    for (int i = 0; i < friction_constraints.size(); i++) {
        CAPTURE(i);
        const FrictionConstraint& constraint = friction_constraints[i];
        const FrictionConstraint& expected_constraint =
            expected_friction_constraints[i];
        if (constraint.closest_point.size() == 1) {
            CHECK(
                constraint.closest_point[0]
                == Catch::Approx(expected_constraint.closest_point[0]));
        } else {
            CHECK(constraint.closest_point.isApprox(
                expected_constraint.closest_point, 1e-12));
        }
        CHECK(constraint.tangent_basis.isApprox(
            expected_constraint.tangent_basis, 1e-12));
        CHECK(
            constraint.normal_force_magnitude
            == Catch::Approx(expected_constraint.normal_force_magnitude));
        CHECK(constraint.mu == Catch::Approx(expected_constraint.mu));

        CHECK(
            constraint.vertex_ids(E, F)
            == expected_constraint.vertex_ids(E, F));
    }

    const Eigen::MatrixXd velocity = V_end - V_start;

    double potential =
        friction_constraints.compute_potential(mesh, velocity, epsv_times_h);

    CHECK(potential == Catch::Approx(expected_potential));

    Eigen::VectorXd grad = friction_constraints.compute_potential_gradient(
        mesh, velocity, epsv_times_h);

    CHECK(grad.isApprox(expected_grad));

    Eigen::SparseMatrix<double> hess =
        friction_constraints.compute_potential_hessian(
            mesh, velocity, epsv_times_h);

    CHECK(hess.isApprox(expected_hess));
}
