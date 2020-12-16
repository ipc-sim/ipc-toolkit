#include <catch2/catch.hpp>

#include <array>
#include <vector>

#include <finitediff.hpp>
#include <igl/edges.h>
#include <igl/readDMAT.h>

#include <ipc/friction/friction.hpp>
#include <ipc/utils/logger.hpp>

#include "../test_utils.hpp"

using namespace ipc;

TEST_CASE("Test friction gradient and hessian", "[friction][gradient][hessian]")
{
    double mu = GENERATE(range(0.0, 1.0, 0.1));
    double epsv_times_h = pow(10, GENERATE(range(-6, 0)));
    double dhat = pow(10, GENERATE(range(-4, 0)));
    double barrier_stiffness = pow(10, GENERATE(range(0, 2)));
    CAPTURE(mu, epsv_times_h, dhat, barrier_stiffness);

    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;
    Constraints contact_constraint_set;
    SECTION("point-triangle")
    {
        V0.resize(4, 3);
        double d = GENERATE_COPY(range(0.0, 2 * dhat, 2 * dhat / 10.0));
        V0.row(0) << 0, d, 0;   // point at t=0
        V0.row(1) << -1, 0, 1;  // triangle vertex 0 at t=0
        V0.row(2) << 2, 0, 0;   // triangle vertex 1 at t=0
        V0.row(3) << -1, 0, -1; // triangle vertex 2 at t=0

        V1 = V0;
        double dy = GENERATE(-1, 1, 1e-1);
        V1.row(0) << 1, d + dy, 0; // point at t=1

        F.resize(1, 3);
        F << 1, 2, 3;
        igl::edges(F, E);
        REQUIRE(E.rows() == 3);

        contact_constraint_set.fv_constraints.emplace_back(0, 0);
    }
    SECTION("edge-edge")
    {
        V0.resize(4, 3);
        double d = GENERATE_COPY(range(0.0, 2 * dhat, 2 * dhat / 10.0));
        V0.row(0) << -1, d, 0; // edge a vertex 0 at t=0
        V0.row(1) << 1, d, 0;  // edge a vertex 1 at t=0
        V0.row(2) << 0, 0, -1; // edge b vertex 0 at t=0
        V0.row(3) << 0, 0, 1;  // edge b vertex 1 at t=0

        V1 = V0;
        double dy = GENERATE(-1, 1, 1e-1);
        V0.row(0) << 0.5, d, 0; // edge a vertex 0 at t=1
        V0.row(1) << 2.5, d, 0; // edge a vertex 1 at t=1

        E.resize(2, 2);
        E.row(0) << 0, 1;
        E.row(1) << 2, 3;

        contact_constraint_set.ee_constraints.emplace_back(0, 1, 0.0);
    }

    FrictionConstraints friction_constraint_set;
    construct_friction_constraint_set(
        V0, E, F, contact_constraint_set, dhat, barrier_stiffness, mu,
        friction_constraint_set);

    Eigen::VectorXd grad = compute_friction_potential_gradient(
        V0, V1, E, F, friction_constraint_set, epsv_times_h);

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        return compute_friction_potential(
            V0, unflatten(x, V1.cols()), E, F, friction_constraint_set,
            epsv_times_h);
    };
    Eigen::VectorXd fgrad;
    fd::finite_gradient(flatten(V1), f, fgrad);
    CHECK(fd::compare_gradient(grad, fgrad));

    Eigen::MatrixXd hess = compute_friction_potential_hessian(
        V0, V1, E, F, friction_constraint_set, epsv_times_h);
    Eigen::MatrixXd fhess;
    fd::finite_hessian(flatten(V1), f, fhess);
    CHECK(fd::compare_hessian(hess, fhess, 1e-3));
}

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
                assert(-mmcvid[3] >= 1);
                constraints.vv_constraints.back().multiplicity = -mmcvid[3];
                normal_force_magnitudes[i] /= -mmcvid[3];
                constraint = &(constraints.vv_constraints.back());

            } else if (mmcvid[3] < 0) { // Is EV?
                edges.emplace_back(mmcvid[1], mmcvid[2]);
                constraints.ev_constraints.emplace_back(
                    edges.size() - 1, -mmcvid[0] - 1);
                constraints.ev_constraints.back().multiplicity = -mmcvid[3];
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
    const std::string& filename_root,
    Eigen::MatrixXd& V0,
    Eigen::MatrixXd& V1,
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
    Eigen::MatrixXd& hess)
{
    Eigen::MatrixXi E1, F1;
    if (!load_mesh(fmt::format("{}_V0.obj", filename_root), V0, E1, F1)) {
        return false;
    }

    if (!load_mesh(fmt::format("{}_V1.obj", filename_root), V1, E1, F1)) {
        return false;
    }

    // MMVCIDs
    Eigen::MatrixXi mmcvids;
    Eigen::VectorXd lambda;
    Eigen::MatrixXd coords, bases;
    bool success = igl::readDMAT(
        TEST_DATA_DIR + fmt::format("{}_mmcvids.dmat", filename_root), mmcvids);
    success |= igl::readDMAT(
        TEST_DATA_DIR + fmt::format("{}_lambda.dmat", filename_root), lambda);
    success |= igl::readDMAT(
        TEST_DATA_DIR + fmt::format("{}_coords.dmat", filename_root), coords);
    success |= igl::readDMAT(
        TEST_DATA_DIR + fmt::format("{}_bases.dmat", filename_root), bases);
    if (!success) {
        return false;
    }

    mmcvids_to_friction_constraints(
        E, F, mmcvids, lambda, coords, bases, friction_constraint_set);
    mmcvids_to_constraints(E, F, mmcvids, constraint_set);

    Eigen::VectorXd params;
    if (!igl::readDMAT(
            TEST_DATA_DIR + fmt::format("{}_params.dmat", filename_root),
            params)
        || params.size() != 4) {
        return false;
    }
    dhat = sqrt(params[0]);         // dHat
    barrier_stiffness = params[1];  // mu_IP
    epsv_times_h = sqrt(params[2]); // fricDHat
    mu = params[3];                 // selfFric

    for (int i = 0; i < friction_constraint_set.size(); i++) {
        friction_constraint_set[i].mu = mu;
    }

    Eigen::VectorXd energy_vec;
    if (!igl::readDMAT(
            TEST_DATA_DIR + fmt::format("{}_energy.dmat", filename_root),
            energy_vec)) {
        return false;
    }
    potential = energy_vec[0];

    if (!igl::readDMAT(
            TEST_DATA_DIR + fmt::format("{}_grad.dmat", filename_root), grad)) {
        return false;
    }

    if (!igl::readDMAT(
            TEST_DATA_DIR + fmt::format("{}_hess.dmat", filename_root), hess)) {
        return false;
    }

    return true;
}

TEST_CASE("Compare IPC friction gradient", "[friction][gradient][hessian]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;
    Constraints contact_constraint_set;
    FrictionConstraints expected_friction_constraint_set;
    double dhat, barrier_stiffness, epsv_times_h, mu;
    double expected_potential;
    Eigen::VectorXd expected_grad;
    Eigen::MatrixXd expected_hess;

    int frame = GENERATE(range(0, 373));
    CAPTURE(frame);

    bool success = read_ipc_friction_data(
        fmt::format("friction/cube_cube_{:d}", frame), V0, V1, E, F,
        contact_constraint_set, expected_friction_constraint_set, dhat,
        barrier_stiffness, epsv_times_h, mu, expected_potential, expected_grad,
        expected_hess);
    REQUIRE(success);

    FrictionConstraints friction_constraint_set;
    construct_friction_constraint_set(
        V0, E, F, contact_constraint_set, dhat, barrier_stiffness, mu,
        friction_constraint_set);

    REQUIRE(friction_constraint_set.size() == contact_constraint_set.size());
    REQUIRE(
        friction_constraint_set.size()
        == expected_friction_constraint_set.size());

    for (int i = 0; i < friction_constraint_set.size(); i++) {
        const FrictionConstraint& constraint = friction_constraint_set[i];
        const FrictionConstraint& expected_constraint =
            expected_friction_constraint_set[i];
        if (constraint.closest_point.size() == 1) {
            CHECK(
                constraint.closest_point[0]
                == Approx(expected_constraint.closest_point[0]));
        } else {
            CHECK(fd::compare_gradient(
                constraint.closest_point, expected_constraint.closest_point,
                1e-12));
        }
        CHECK(fd::compare_jacobian(
            constraint.tangent_basis, expected_constraint.tangent_basis,
            1e-12));
        CHECK(
            constraint.normal_force_magnitude
            == Approx(expected_constraint.normal_force_magnitude));
        CHECK(constraint.mu == Approx(expected_constraint.mu));
    }

    double potential = compute_friction_potential(
        V0, V1, E, F, friction_constraint_set, epsv_times_h);

    CHECK(potential == Approx(expected_potential));

    Eigen::VectorXd grad = compute_friction_potential_gradient(
        V0, V1, E, F, friction_constraint_set, epsv_times_h);

    CHECK(fd::compare_gradient(expected_grad, grad, 1e-8));

    Eigen::SparseMatrix<double> hess = compute_friction_potential_hessian(
        V0, V1, E, F, friction_constraint_set, epsv_times_h);

    CHECK(fd::compare_hessian(expected_hess, Eigen::MatrixXd(hess), 1e-6));
}
