#include <catch2/catch.hpp>

#include <array>
#include <vector>

#include <finitediff.hpp>
#include <igl/edges.h>

#include <ipc/ipc.hpp>
#include <ipc/friction/friction.hpp>
#include <ipc/utils/logger.hpp>
#include <ipc/config.hpp>

#include "friction_data_generator.hpp"
#include "../test_utils.hpp"

#include <unsupported/Eigen/SparseExtra>

using namespace ipc;

TEST_CASE("Test friction gradient and hessian", "[friction][gradient][hessian]")
{
    FrictionData data = friction_data_generator();
    const auto& [V0, V1, E, F, contact_constraint_set, mu, epsv_times_h, dhat, barrier_stiffness] =
        data;

    CollisionMesh mesh(V0, E, F);

    FrictionConstraints friction_constraint_set;
    construct_friction_constraint_set(
        mesh, V0, contact_constraint_set, dhat, barrier_stiffness, mu,
        friction_constraint_set);

    Eigen::VectorXd grad = compute_friction_potential_gradient(
        mesh, V0, V1, friction_constraint_set, epsv_times_h);

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        return compute_friction_potential(
            mesh, data.V0, fd::unflatten(x, data.V1.cols()),
            friction_constraint_set, data.epsv_times_h);
    };
    Eigen::VectorXd fgrad;
    fd::finite_gradient(fd::flatten(V1), f, fgrad);
    CHECK(fd::compare_gradient(grad, fgrad));

    Eigen::MatrixXd hess = compute_friction_potential_hessian(
        mesh, V0, V1, friction_constraint_set, epsv_times_h);
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
    Constraints& constraint_set,
    FrictionConstraints& friction_constraint_set,
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
        E, F, mmcvids, lambda, coords, bases, friction_constraint_set);
    mmcvids_to_constraints(E, F, mmcvids, constraint_set);
    for (int i = 0; i < friction_constraint_set.size(); i++) {
        friction_constraint_set[i].mu = mu;
    }

    return true;
}

TEST_CASE(
    "Compare IPC friction derivatives", "[friction][gradient][hessian][data]")
{
    Eigen::MatrixXd V_start, V_lagged, V_end;
    Eigen::MatrixXi E, F;
    Constraints contact_constraint_set;
    FrictionConstraints expected_friction_constraint_set;
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
        V_start, V_lagged, V_end, E, F, contact_constraint_set,
        expected_friction_constraint_set, dhat, barrier_stiffness, epsv_times_h,
        mu, expected_potential, expected_grad, expected_hess);
    REQUIRE(success);

    Eigen::MatrixXi face_edges;
    igl::edges(F, face_edges);
    E.conservativeResize(E.rows() + face_edges.rows(), E.cols());
    E.bottomRows(face_edges.rows()) = face_edges;
    CollisionMesh mesh(V_start, E, F);

    FrictionConstraints friction_constraint_set;
    construct_friction_constraint_set(
        mesh, V_lagged, contact_constraint_set, dhat, barrier_stiffness, mu,
        friction_constraint_set);

    REQUIRE(friction_constraint_set.size() == contact_constraint_set.size());
    REQUIRE(
        friction_constraint_set.size()
        == expected_friction_constraint_set.size());

    REQUIRE(V_start.size() == V_lagged.size());
    REQUIRE(V_start.size() == V_end.size());
    REQUIRE(V_start.size() == expected_grad.size());

    for (int i = 0; i < friction_constraint_set.size(); i++) {
        CAPTURE(i);
        const FrictionConstraint& constraint = friction_constraint_set[i];
        const FrictionConstraint& expected_constraint =
            expected_friction_constraint_set[i];
        if (constraint.closest_point.size() == 1) {
            CHECK(
                constraint.closest_point[0]
                == Approx(expected_constraint.closest_point[0]));
        } else {
            CHECK(constraint.closest_point.isApprox(
                expected_constraint.closest_point, 1e-12));
        }
        CHECK(constraint.tangent_basis.isApprox(
            expected_constraint.tangent_basis, 1e-12));
        CHECK(
            constraint.normal_force_magnitude
            == Approx(expected_constraint.normal_force_magnitude));
        CHECK(constraint.mu == Approx(expected_constraint.mu));

        CHECK(
            constraint.vertex_indices(E, F)
            == expected_constraint.vertex_indices(E, F));
    }

    double potential = compute_friction_potential(
        mesh, V_start, V_end, friction_constraint_set, epsv_times_h);

    CHECK(potential == Approx(expected_potential));

    Eigen::VectorXd grad = compute_friction_potential_gradient(
        mesh, V_start, V_end, friction_constraint_set, epsv_times_h);

    CHECK(grad.isApprox(expected_grad));

    Eigen::SparseMatrix<double> hess = compute_friction_potential_hessian(
        mesh, V_start, V_end, friction_constraint_set, epsv_times_h);

    CHECK(hess.isApprox(expected_hess));
}
