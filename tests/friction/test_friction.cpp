#include <catch2/catch.hpp>

#include <array>
#include <vector>

#include <finitediff.hpp>
#include <igl/edges.h>
#include <nlohmann/json.hpp>

#include <ipc/ipc.hpp>
#include <ipc/friction/friction.hpp>
#include <ipc/utils/logger.hpp>

#include "friction_data_generator.hpp"
#include "../test_utils.hpp"

#include <unsupported/Eigen/SparseExtra>

using namespace ipc;

TEST_CASE("Test friction gradient and hessian", "[friction][gradient][hessian]")
{
    const FrictionData& data = GENERATE(FrictionDataGenerator::create());
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

TEST_CASE("Test friction force jacobian", "[friction][force-jacobian]")
{
    const FrictionData& data = GENERATE(FrictionDataGenerator::create());
    const auto& [V0, V1, E, F, contact_constraint_set, mu, epsv_times_h, dhat, barrier_stiffness] =
        data;

    CollisionMesh mesh(V0, E, F);

    double distance_t0 =
        compute_minimum_distance(mesh, V0, contact_constraint_set);
    double distance_t1 =
        compute_minimum_distance(mesh, V1, contact_constraint_set);
    // CHECK((distance_t0 < dhat || distance_t1 < dhat));

    if (distance_t0 == 0 || distance_t1 == 0) {
        return;
    }

    CAPTURE(
        data.mu, data.epsv_times_h, data.dhat, data.barrier_stiffness,
        contact_constraint_set.size());

    Eigen::MatrixXd X, Ut, U;
    SECTION("X = V0") { X = V0; }
    SECTION("X = V0 - (V1 - V0)") { X = V0 - (V1 - V0); }
    Ut = V0 - X;
    U = V1 - X;

    FrictionConstraints friction_constraint_set;
    construct_friction_constraint_set(
        mesh, X + Ut, contact_constraint_set, dhat, barrier_stiffness, mu,
        friction_constraint_set);

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_X = compute_friction_force_jacobian(
        mesh, X, Ut, U, friction_constraint_set, dhat, barrier_stiffness,
        epsv_times_h, FrictionConstraint::DiffWRT::X);

    auto F_X = [&](const Eigen::VectorXd& x) {
        construct_friction_constraint_set(
            mesh, fd::unflatten(x, X.cols()) + Ut, data.constraints, data.dhat,
            data.barrier_stiffness, data.mu, friction_constraint_set);
        return compute_friction_force(
            mesh, fd::unflatten(x, X.cols()), Ut, U, friction_constraint_set,
            data.dhat, data.barrier_stiffness, data.epsv_times_h);
    };
    Eigen::MatrixXd fd_JF_wrt_X;
    fd::finite_jacobian(fd::flatten(X), F_X, fd_JF_wrt_X);
    CHECK(fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X));

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_Ut = compute_friction_force_jacobian(
        mesh, X, Ut, U, friction_constraint_set, dhat, barrier_stiffness,
        epsv_times_h, FrictionConstraint::DiffWRT::Ut);

    auto F_Ut = [&](const Eigen::VectorXd& ut) {
        construct_friction_constraint_set(
            mesh, X + fd::unflatten(ut, Ut.cols()), data.constraints, data.dhat,
            data.barrier_stiffness, data.mu, friction_constraint_set);
        return compute_friction_force(
            mesh, X, fd::unflatten(ut, Ut.cols()), U, friction_constraint_set,
            data.dhat, data.barrier_stiffness, data.epsv_times_h);
    };
    Eigen::MatrixXd fd_JF_wrt_Ut;
    fd::finite_jacobian(fd::flatten(Ut), F_Ut, fd_JF_wrt_Ut);
    CHECK(fd::compare_jacobian(JF_wrt_Ut, fd_JF_wrt_Ut));

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_U = compute_friction_force_jacobian(
        mesh, X, Ut, U, friction_constraint_set, dhat, barrier_stiffness,
        epsv_times_h, FrictionConstraint::DiffWRT::U);

    auto F_U = [&](const Eigen::VectorXd& u) {
        return compute_friction_force(
            mesh, X, Ut, fd::unflatten(u, U.cols()), friction_constraint_set,
            data.dhat, data.barrier_stiffness, data.epsv_times_h);
    };
    Eigen::MatrixXd fd_JF_wrt_U;
    fd::finite_jacobian(fd::flatten(U), F_U, fd_JF_wrt_U);
    CHECK(fd::compare_jacobian(JF_wrt_U, fd_JF_wrt_U));

    ///////////////////////////////////////////////////////////////////////////

    Eigen::VectorXd force = compute_friction_force(
        mesh, X, Ut, U, friction_constraint_set, dhat, barrier_stiffness,
        epsv_times_h);
    Eigen::VectorXd grad_D = compute_friction_potential_gradient(
        mesh, X + Ut, X + U, friction_constraint_set, epsv_times_h);
    CHECK(fd::compare_gradient(-force, grad_D));

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd jac_force = compute_friction_force_jacobian(
        mesh, X, Ut, U, friction_constraint_set, dhat, barrier_stiffness,
        epsv_times_h, FrictionConstraint::DiffWRT::U);
    Eigen::MatrixXd hess_D = compute_friction_potential_hessian(
        mesh, X + Ut, X + U, friction_constraint_set, epsv_times_h, false);
    CHECK(fd::compare_jacobian(-jac_force, hess_D));
}

inline Eigen::MatrixXd loadMarketXd(const std::string& f)
{
    Eigen::SparseMatrix<double> tmp;
    REQUIRE(Eigen::loadMarket(tmp, f));
    return Eigen::MatrixXd(tmp);
}

inline Eigen::MatrixXi loadMarketXi(const std::string& f)
{
    Eigen::SparseMatrix<int> tmp;
    REQUIRE(Eigen::loadMarket(tmp, f));
    return Eigen::MatrixXi(tmp);
}

TEST_CASE(
    "Test friction force jacobian on real data",
    "[friction][force-jacobian][thisone]")
{
    std::string scene;
    bool is_2D = true;
    SECTION("point-plane") { scene = "point-plane"; }
    // SECTION("square-incline") { scene = "square-incline"; }

    Eigen::MatrixXd X, Ut, U;
    Eigen::MatrixXi E, F;
    {
        X = loadMarketXd(fmt::format(
            "{}friction-force-jacobian/{}/X.mat", TEST_DATA_DIR, scene));
        Ut = loadMarketXd(fmt::format(
            "{}friction-force-jacobian/{}/Ut.mat", TEST_DATA_DIR, scene));
        Ut = fd::unflatten(Ut, X.cols());
        REQUIRE(X.rows() == Ut.rows());
        REQUIRE(X.cols() == Ut.cols());
        U = loadMarketXd(fmt::format(
            "{}friction-force-jacobian/{}/U.mat", TEST_DATA_DIR, scene));
        U = fd::unflatten(U, X.cols());
        REQUIRE(X.rows() == U.rows());
        REQUIRE(X.cols() == U.cols());
        if (is_2D) {
            E = loadMarketXi(fmt::format(
                "{}friction-force-jacobian/{}/F.mat", TEST_DATA_DIR, scene));
        } else {
            F = loadMarketXi(fmt::format(
                "{}friction-force-jacobian/{}/F.mat", TEST_DATA_DIR, scene));
            igl::edges(F, E);
        }
    }

    std::vector<bool> is_on_surface =
        CollisionMesh::construct_is_on_surface(X.rows(), E);
    CollisionMesh mesh(is_on_surface, X, E, F);

    X = mesh.vertices(X);
    Ut = mesh.vertices(Ut);
    U = mesh.vertices(U);

    E = mesh.edges();
    F = mesh.faces();

    double mu = 0.5;
    double dhat = 0.1;
    double kappa = 141;
    double epsv_dt = 5e-6;

    Constraints contact_constraint_set;
    construct_constraint_set(mesh, X + Ut, dhat, contact_constraint_set);

    CHECK(compute_minimum_distance(mesh, X + Ut, contact_constraint_set) != 0);
    CHECK(compute_minimum_distance(mesh, X + U, contact_constraint_set) != 0);

    FrictionConstraints friction_constraint_set;
    construct_friction_constraint_set(
        mesh, X + Ut, contact_constraint_set, dhat, kappa, mu,
        friction_constraint_set);

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_X = compute_friction_force_jacobian(
        mesh, X, Ut, U, friction_constraint_set, dhat, kappa, epsv_dt,
        FrictionConstraint::DiffWRT::X);

    auto F_X = [&](const Eigen::VectorXd& x) {
        CollisionMesh mesh(fd::unflatten(x, X.cols()), E, F);
        Constraints contact_constraint_set;
        construct_constraint_set(
            mesh, fd::unflatten(x, X.cols()) + Ut, dhat,
            contact_constraint_set);
        FrictionConstraints friction_constraint_set;
        construct_friction_constraint_set(
            mesh, fd::unflatten(x, X.cols()) + Ut, contact_constraint_set, dhat,
            kappa, mu, friction_constraint_set);
        return compute_friction_force(
            mesh, fd::unflatten(x, X.cols()), Ut, U, friction_constraint_set,
            dhat, kappa, epsv_dt);
    };
    Eigen::MatrixXd fd_JF_wrt_X;
    fd::finite_jacobian(fd::flatten(X), F_X, fd_JF_wrt_X);
    CHECK(fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X));
    // std::cout << JF_wrt_X << "\n\n" << fd_JF_wrt_X << std::endl;

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_Ut = compute_friction_force_jacobian(
        mesh, X, Ut, U, friction_constraint_set, dhat, kappa, epsv_dt,
        FrictionConstraint::DiffWRT::Ut);

    auto F_Ut = [&](const Eigen::VectorXd& ut) {
        Constraints contact_constraint_set;
        construct_constraint_set(
            mesh, X + fd::unflatten(ut, Ut.cols()), dhat,
            contact_constraint_set);
        FrictionConstraints friction_constraint_set;
        construct_friction_constraint_set(
            mesh, X + fd::unflatten(ut, Ut.cols()), contact_constraint_set,
            dhat, kappa, mu, friction_constraint_set);
        return compute_friction_force(
            mesh, X, fd::unflatten(ut, Ut.cols()), U, friction_constraint_set,
            dhat, kappa, epsv_dt);
    };
    Eigen::MatrixXd fd_JF_wrt_Ut;
    fd::finite_jacobian(fd::flatten(Ut), F_Ut, fd_JF_wrt_Ut);
    CAPTURE(scene);
    CHECK(fd::compare_jacobian(JF_wrt_Ut, fd_JF_wrt_Ut));
    // std::cout << JF_wrt_Ut << "\n\n" << fd_JF_wrt_Ut << "\n" << std::endl;

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_U = compute_friction_force_jacobian(
        mesh, X, Ut, U, friction_constraint_set, dhat, kappa, epsv_dt,
        FrictionConstraint::DiffWRT::U);

    auto F_U = [&](const Eigen::VectorXd& u) {
        return compute_friction_force(
            mesh, X, Ut, fd::unflatten(u, U.cols()), friction_constraint_set,
            dhat, kappa, epsv_dt);
    };
    Eigen::MatrixXd fd_JF_wrt_U;
    fd::finite_jacobian(fd::flatten(U), F_U, fd_JF_wrt_U);

    CHECK(fd::compare_jacobian(JF_wrt_U, fd_JF_wrt_U));
    // std::cout << JF_wrt_U << "\n\n" << fd_JF_wrt_U << "\n" << std::endl;

    ///////////////////////////////////////////////////////////////////////////

    Eigen::VectorXd force = compute_friction_force(
        mesh, X, Ut, U, friction_constraint_set, dhat, kappa, epsv_dt);
    Eigen::VectorXd grad_D = compute_friction_potential_gradient(
        mesh, X + Ut, X + U, friction_constraint_set, epsv_dt);
    CHECK(fd::compare_gradient(-force, grad_D));

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd jac_force = compute_friction_force_jacobian(
        mesh, X, Ut, U, friction_constraint_set, dhat, kappa, epsv_dt,
        FrictionConstraint::DiffWRT::U);
    Eigen::MatrixXd hess_D = compute_friction_potential_hessian(
        mesh, X + Ut, X + U, friction_constraint_set, epsv_dt, false);
    CHECK(fd::compare_jacobian(-jac_force, hess_D));
    // std::cout << (-jac_force - hess_D).lpNorm<Eigen::Infinity>() <<
    // std::endl;
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
                constraints.vv_constraints.back().multiplicity() = -mmcvid[3];
                normal_force_magnitudes[i] /= -mmcvid[3];
                constraint = &(constraints.vv_constraints.back());

            } else if (mmcvid[3] < 0) { // Is EV?
                edges.emplace_back(mmcvid[1], mmcvid[2]);
                constraints.ev_constraints.emplace_back(
                    edges.size() - 1, -mmcvid[0] - 1);
                CHECK(-mmcvid[3] >= 1);
                constraints.ev_constraints.back().multiplicity() = -mmcvid[3];
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

template <typename T>
inline void from_json(const nlohmann::json& json, VectorX<T>& vec)
{
    vec =
        Eigen::Map<VectorX<T>>(json.get<std::vector<T>>().data(), json.size());
}

template <typename T>
void from_json(const nlohmann::json& json, MatrixX<T>& mat)
{
    typedef std::vector<std::vector<T>> L;
    L list = json.get<L>();

    size_t num_rows = list.size();
    if (num_rows == 0) {
        return;
    }
    size_t num_cols = list[0].size();
    mat.resize(num_rows, num_cols);

    for (size_t i = 0; i < num_rows; ++i) {
        assert(num_cols == list[i].size());
        mat.row(i) = Eigen::Map<RowVectorX<T>>(list[i].data(), num_cols);
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
        IPC_LOG(error("Unable to open IPC friction data file: {}", filename));
        return false;
    }

    if (data.is_discarded()) {
        IPC_LOG(error("IPC friction data JSON is invalid: {}", filename));
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
