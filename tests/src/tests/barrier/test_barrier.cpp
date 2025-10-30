#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>

#include <ipc/barrier/barrier.hpp>
#include <ipc/geometry/normal.hpp>
#include <ipc/smooth_contact/primitives/point3.hpp>
#include <ipc/utils/autodiff_types.hpp>
#include <ipc/math/math.hpp>

#include <finitediff.hpp>
#include <igl/edges.h>

#include <memory>

using namespace ipc;

/// @warning This implementation will not work with dmin > 0
class PhysicalBarrier : public NormalizedClampedLogBarrier {
public:
    PhysicalBarrier(const bool _use_dist_sqr) : use_dist_sqr(_use_dist_sqr) { }

    double operator()(const double d, const double dhat) const override
    {
        return (use_dist_sqr ? sqrt(dhat) : dhat)
            * NormalizedClampedLogBarrier::operator()(d, dhat);
    }

    double first_derivative(const double d, const double dhat) const override
    {
        return (use_dist_sqr ? sqrt(dhat) : dhat)
            * NormalizedClampedLogBarrier::first_derivative(d, dhat);
    }

    double second_derivative(const double d, const double dhat) const override
    {
        return (use_dist_sqr ? sqrt(dhat) : dhat)
            * NormalizedClampedLogBarrier::second_derivative(d, dhat);
    }

private:
    bool use_dist_sqr;
};

TEST_CASE("Spline derivatives", "[deriv]")
{
    const int n_samples = 100;
    ScalarBase::setVariableCount(1);
    using T = ipc::ADHessian<1>;
    for (int i = 0; i < n_samples; i++) {
        const double x = 2. * i / static_cast<double>(n_samples) - 1;
        double deriv = ipc::Math<double>::cubic_spline_grad(x);
        double hess = ipc::Math<double>::cubic_spline_hess(x);
        T x_ad = T(x, 0);
        T y_ad = ipc::Math<T>::cubic_spline(x_ad);
        double deriv_ad = y_ad.grad(0);
        double hess_ad = y_ad.Hess(0);

        CHECK(abs(deriv_ad - deriv) < 1e-14 * std::max(1., abs(deriv)));
        CHECK(abs(hess_ad - hess) < 1e-12 * std::max(1., abs(hess)));
    }
}

TEST_CASE("Heaviside derivatives", "[deriv]")
{
    const int n_samples = 100;
    ScalarBase::setVariableCount(1);
    using T = ipc::ADHessian<1>;
    for (int i = 0; i < n_samples; i++) {
        const double x = i / static_cast<double>(n_samples) - 1;
        double deriv = ipc::Math<double>::smooth_heaviside_grad(x, 1., 0.);
        double hess = ipc::Math<double>::smooth_heaviside_hess(x, 1., 0.);
        T x_ad = T(x, 0);
        T y_ad = ipc::Math<T>::smooth_heaviside(x_ad, 1., 0.);
        double deriv_ad = y_ad.grad(0);
        double hess_ad = y_ad.Hess(0);

        CHECK(abs(deriv_ad - deriv) < 1e-14 * std::max(1., abs(deriv)));
        CHECK(abs(hess_ad - hess) < 1e-12 * std::max(1., abs(hess)));
    }
}

TEST_CASE("Inv barrier derivatives", "[deriv]")
{
    const int n_samples = 100;
    ScalarBase::setVariableCount(1);
    using T = ipc::ADHessian<1>;
    const int r = 1;
    const double dhat = 0.13;
    for (int i = 1; i <= n_samples; i++) {
        const double x = i / static_cast<double>(n_samples);
        double deriv = ipc::Math<double>::inv_barrier_grad(x / dhat, r) / dhat;
        double hess =
            ipc::Math<double>::inv_barrier_hess(x / dhat, r) / dhat / dhat;
        T x_ad = T(x, 0);
        T y_ad = ipc::Math<T>::inv_barrier(x_ad / dhat, r);
        double deriv_ad = y_ad.grad(0);
        double hess_ad = y_ad.Hess(0);

        CHECK(abs(deriv_ad - deriv) < 1e-14 * std::max(1., abs(deriv)));
        CHECK(abs(hess_ad - hess) < 1e-12 * std::max(1., abs(hess)));
    }

    ScalarBase::setVariableCount(3);
    using T3 = ipc::ADHessian<3>;

    for (int i = 1; i <= n_samples; i++) {
        Eigen::Vector3d x = Eigen::Vector3d::Random() * dhat / 3.;
        double deriv =
            ipc::Math<double>::inv_barrier_grad(x.norm() / dhat, r) / dhat;
        double hess = ipc::Math<double>::inv_barrier_hess(x.norm() / dhat, r)
            / dhat / dhat;
        auto x_ad = ipc::slice_positions<T3, 3, 1>(x);
        T3 y_ad = ipc::Math<T3>::inv_barrier(x_ad.norm() / dhat, r);
        Eigen::Vector3d deriv_ad = y_ad.grad;
        Eigen::Matrix3d hess_ad = y_ad.Hess;

        Eigen::Vector3d xn = x / x.norm();
        Eigen::Vector3d deriv_analytic = deriv * xn;
        Eigen::Matrix3d hess_analytic =
            (deriv / x.norm()) * Eigen::Matrix3d::Identity()
            + (hess - deriv / x.norm()) * xn * xn.transpose();

        CHECK(
            (deriv_ad - deriv_analytic).norm()
            < 1e-14 * std::max(1., deriv_ad.norm()));
        CHECK(
            (hess_ad - hess_analytic).norm()
            < 1e-12 * std::max(1., hess_ad.norm()));
    }
}

TEST_CASE("Normalize vector derivatives", "[deriv]")
{
    const int n_samples = 1000;
    ScalarBase::setVariableCount(3);
    using T = ipc::ADHessian<3>;
    for (int i = 1; i <= n_samples; i++) {
        Eigen::Vector3d x = Eigen::Vector3d::Random();
        ipc::Vector<T, 3> x_ad = ipc::slice_positions<T, 3, 1>(x);
        ipc::Vector<T, 3> y_ad = x_ad / x_ad.norm();

        const auto [y, grad, hess] =
            ipc::normalization_and_jacobian_and_hessian(x);

        double err_grad = 0, err_hess = 0;
        for (int j = 0; j < 3; j++) {
            err_grad += (grad.col(j) - y_ad(j).grad).norm();
            err_hess += (hess[j] - y_ad(j).Hess).norm();
        }

        CHECK((grad - grad.transpose()).norm() <= 1e-12);
        CHECK((hess[0] - hess[0].transpose()).norm() <= 1e-12);
        CHECK(err_grad <= 1e-12);
        CHECK(err_hess <= 1e-12);
    }
}

TEST_CASE("line-line closest direction derivatives", "[deriv]")
{
    const int n_samples = 100;
    ScalarBase::setVariableCount(12);
    using T = ipc::ADHessian<12>;
    for (int i = 1; i <= n_samples; i++) {
        ipc::Vector6d ea = ipc::Vector6d::Random();
        ipc::Vector<T, 6> eaT = ipc::slice_positions<T, 6, 1>(ea);
        ipc::Vector6d eb = ipc::Vector6d::Random();
        ipc::Vector<T, 6> ebT = ipc::slice_positions<T, 6, 1>(eb, 6);

        ipc::Vector<T, 3> dT = ipc::line_line_closest_point_direction<T>(
            eaT.head<3>(), eaT.tail<3>(), ebT.head<3>(), ebT.tail<3>());

        const auto [d1, grad1] =
            ipc::line_line_closest_point_direction_gradient(
                ea.head<3>(), ea.tail<3>(), eb.head<3>(), eb.tail<3>());

        const auto [d2, grad2, hess2] =
            ipc::line_line_closest_point_direction_hessian(
                ea.head<3>(), ea.tail<3>(), eb.head<3>(), eb.tail<3>());

        double err_grad = (grad1 - grad2).norm();
        double err_hess = 0;
        for (int j = 0; j < 3; j++) {
            err_grad += (grad2.row(j).transpose() - dT(j).grad).norm();
            err_hess += (hess2[j] - dT(j).Hess).norm();
        }

        CHECK(err_grad <= 1e-12);
        CHECK(err_hess <= 1e-10);
    }
}

TEST_CASE("opposite_direction_penalty derivatives", "[deriv]")
{
    const int n_samples = 1000;
    ScalarBase::setVariableCount(6);
    using T = ipc::ADHessian<6>;
    const double alpha = 1;
    const double beta = 1;
    for (int i = 1; i <= n_samples; i++) {
        ipc::Vector6d x = ipc::Vector6d::Random();
        ipc::Vector<T, 6> x_ad = ipc::slice_positions<T, 6, 1>(x);
        T y_ad = ipc::Math<T>::smooth_heaviside(
            x_ad.tail(3).dot(x_ad.head(3)) / x_ad.head(3).norm(), alpha, beta);

        const auto [y, grad, hess] = ipc::opposite_direction_penalty_hess(
            x.head(3), x.tail(3), alpha, beta);

        CHECK((grad - y_ad.grad).norm() <= 1e-12);
        CHECK((hess - y_ad.Hess).norm() <= 1e-12);
    }
}

TEST_CASE("negative_orientation_penalty derivatives", "[deriv]")
{
    const int n_samples = 1000;
    const double alpha = 1;
    const double beta = 1;
    for (int i = 1; i <= n_samples; i++) {
        ScalarBase::setVariableCount(9);
        using T = ipc::ADHessian<9>;
        ipc::Vector9d x = ipc::Vector9d::Random();
        ipc::Vector<T, 9> x_ad = ipc::slice_positions<T, 9, 1>(x);
        ipc::Vector<T, 3> t = x_ad.head<3>().cross(x_ad.segment<3>(3));
        T y_ad = ipc::Math<T>::smooth_heaviside(
            x_ad.tail(3).dot(t) / t.norm(), alpha, beta);

        auto [y, grad, hess] = ipc::negative_orientation_penalty_hess(
            x.head(3), x.segment(3, 3), x.tail(3), alpha, beta);

        // hess.topLeftCorner(6, 6).setZero();
        // auto hess_ad = y_ad.Hess;
        // hess_ad.topLeftCorner(6, 6).setZero();
        // y_ad = T(y_ad.val, y_ad.grad, hess_ad);

        CHECK(abs(y - y_ad.val) <= 1e-12);
        CHECK((grad - y_ad.grad).norm() <= 1e-12);
        CHECK((hess - y_ad.Hess).norm() <= 1e-10);
    }
}

TEST_CASE("point term derivatives", "[deriv]")
{
    ipc::SmoothContactParameters params(1, 1, 1, 0.01, 0, 2);

    Eigen::Matrix<double, -1, 3> vectors(9, 3);
    vectors << -0.696515, -0.173578, -0.696231, 0.50146, -0.0017947, 0.999718,
        0.346908, 0.152939, 1.00049, 0.208062, 0.290611, 1.00075, 0.280725,
        0.498548, 1.00069, 0.499796, 0.49999, 1.00127, 0.501162, 0.498839,
        0.779581, 0.50047, 0.290411, 0.709051, 0.500615, 0.153465, 0.846423;

    Eigen::Matrix<double, -1, 3> V = vectors;
    V.row(0).setZero();

    Eigen::MatrixXi E(8, 2), F(8, 3);
    E << 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8;
    F << 0, 2, 1, 0, 3, 2, 0, 4, 3, 0, 5, 4, 0, 6, 5, 0, 6, 7, 0, 7, 8, 0, 8, 1;
    igl::edges(F, E);
    ipc::CollisionMesh mesh(V, E, F);

    auto point_term =
        std::make_unique<ipc::Point3>(0, mesh, V, -vectors.row(0), params);

    {
        auto [y, y_grad, y_hess] =
            point_term->smooth_point3_term_tangent_hessian(
                vectors.row(0), V.bottomRows(8), params.alpha_t, params.beta_t);

        Eigen::VectorXd fgrad;
        fd::finite_gradient(
            fd::flatten(vectors),
            [&](const Eigen::VectorXd& x) {
                Eigen::MatrixXd V_fd = fd::unflatten(x, 3);
                V_fd.row(0).setZero();
                auto point_term_fd = std::make_unique<ipc::Point3>(
                    0, mesh, V_fd, -x.head<3>(), params);
                return std::get<0>(
                    point_term_fd->smooth_point3_term_tangent_gradient(
                        x.head(3), fd::unflatten(x.tail(x.size() - 3), 3),
                        params.alpha_t, params.beta_t));
            },
            fgrad);

        CHECK((y_grad - fgrad).norm() <= 1e-7 * fgrad.norm());

        Eigen::MatrixXd fhess;
        fd::finite_jacobian(
            fd::flatten(vectors),
            [&](const Eigen::VectorXd& x) {
                Eigen::MatrixXd V_fd = fd::unflatten(x, 3);
                V_fd.row(0).setZero();
                auto point_term_fd = std::make_unique<ipc::Point3>(
                    0, mesh, V_fd, -x.head<3>(), params);
                return std::get<1>(
                    point_term_fd->smooth_point3_term_tangent_gradient(
                        x.head(3), fd::unflatten(x.tail(x.size() - 3), 3),
                        params.alpha_t, params.beta_t));
            },
            fhess);

        CHECK((y_hess - fhess).norm() <= 1e-7 * fhess.norm());
    }

    {
        auto y_grad = point_term->grad(vectors.row(0), fd::flatten(V));
        auto y_hess = point_term->hessian(vectors.row(0), fd::flatten(V));

        Eigen::MatrixXd X(V.rows() + 1, 3);
        X << vectors.row(0), V;

        Eigen::VectorXd fgrad;
        fd::finite_gradient(
            fd::flatten(X),
            [&](const Eigen::VectorXd& x) {
                return point_term->potential(x.head(3), x.tail(x.size() - 3));
            },
            fgrad);

        CHECK((y_grad - fgrad).norm() <= 1e-7 * fgrad.norm());

        Eigen::MatrixXd fhess;
        fd::finite_jacobian(
            fd::flatten(X),
            [&](const Eigen::VectorXd& x) {
                return point_term->grad(x.head(3), x.tail(x.size() - 3));
            },
            fhess);

        CHECK((y_hess - fhess).norm() <= 1e-7 * fhess.norm());
    }
}

TEST_CASE("point term normal derivatives", "[deriv]")
{
    ipc::SmoothContactParameters params(1, 1, 1, 1, 0, 2);

    Eigen::Matrix<double, -1, 3> vectors(9, 3);
    vectors << -0.696515, -0.173578, -0.696231, 0.50146, -0.0017947, 0.999718,
        0.346908, 0.152939, 1.00049, 0.208062, 0.290611, 1.00075, 0.280725,
        0.498548, 1.00069, 0.499796, 0.49999, 1.00127, 0.501162, 0.498839,
        0.779581, 0.50047, 0.290411, 0.709051, 0.500615, 0.153465, 0.846423;

    Eigen::Matrix<double, -1, 3> V = vectors;
    V.row(0).setZero();

    Eigen::MatrixXi E(3, 2), F(1, 3);
    E << 0, 1, 0, 2, 1, 2;
    F << 0, 1, 2;
    igl::edges(F, E);
    ipc::CollisionMesh mesh(V, E, F);

    auto point_term =
        std::make_unique<ipc::Point3>(0, mesh, V, -vectors.row(0), params);

    {
        auto [y, y_grad, y_hess] =
            point_term->smooth_point3_term_normal_hessian(
                vectors.row(0), V.bottomRows(8), params.alpha_n, params.beta_n);

        Eigen::VectorXd fgrad;
        fd::finite_gradient(
            fd::flatten(vectors),
            [&](const Eigen::VectorXd& x) {
                Eigen::MatrixXd V_fd = fd::unflatten(x, 3);
                V_fd.row(0).setZero();
                auto point_term_fd = std::make_unique<ipc::Point3>(
                    0, mesh, V_fd, -x.head<3>(), params);
                return std::get<0>(
                    point_term_fd->smooth_point3_term_normal_hessian(
                        x.head(3), fd::unflatten(x.tail(x.size() - 3), 3),
                        params.alpha_n, params.beta_n));
            },
            fgrad);

        CHECK((y_grad - fgrad).norm() <= 1e-7 * fgrad.norm());

        Eigen::MatrixXd fhess;
        fd::finite_jacobian(
            fd::flatten(vectors),
            [&](const Eigen::VectorXd& x) {
                Eigen::MatrixXd V_fd = fd::unflatten(x, 3);
                V_fd.row(0).setZero();
                auto point_term_fd = std::make_unique<ipc::Point3>(
                    0, mesh, V_fd, -x.head<3>(), params);
                return std::get<1>(
                    point_term_fd->smooth_point3_term_normal_hessian(
                        x.head(3), fd::unflatten(x.tail(x.size() - 3), 3),
                        params.alpha_n, params.beta_n));
            },
            fhess, fd::AccuracyOrder::FOURTH, 1e-6);

        CHECK((y_hess - fhess).norm() <= 1e-7 * fhess.norm());
    }
}

TEST_CASE("Barrier derivatives", "[barrier]")
{
    bool use_dist_sqr = GENERATE(false, true);
    double dhat = GENERATE_COPY(range(use_dist_sqr ? -2 : -5, 0));
    dhat = pow(10, dhat);

    double d =
        GENERATE_COPY(take(10, random(dhat / 2, 0.9 * dhat))); // ∈ [0, d̂]
    Eigen::Matrix<double, 1, 1> d_vec;
    d_vec << d;

    // Check derivatives

    std::unique_ptr<ipc::Barrier> barrier;
    SECTION("Original")
    {
        barrier = std::make_unique<ipc::ClampedLogBarrier>();
    }
    SECTION("Normalized")
    {
        barrier = std::make_unique<ipc::NormalizedClampedLogBarrier>();
    }
    SECTION("Physical")
    {
        barrier = std::make_unique<PhysicalBarrier>(use_dist_sqr);
    }
    SECTION("ClampedLogSq")
    {
        barrier = std::make_unique<ipc::ClampedLogSqBarrier>();
    }
    SECTION("Cubic") { barrier = std::make_unique<ipc::CubicBarrier>(); }
    SECTION("TwoStage") { barrier = std::make_unique<ipc::TwoStageBarrier>(); }

    if (use_dist_sqr) {
        d_vec *= d;
        d *= d;
        dhat *= dhat;
    }

    Eigen::VectorXd fgrad(1);
    fd::finite_gradient(
        d_vec,
        [&](const Eigen::VectorXd& _d) { return (*barrier)(_d[0], dhat); },
        fgrad);

    Eigen::VectorXd grad(1);
    grad << barrier->first_derivative(d, dhat);

    CAPTURE(dhat, d, fgrad(0), grad(0), use_dist_sqr);
    CHECK(fd::compare_gradient(fgrad, grad));

    // Check second_derivative

    fd::finite_gradient(
        d_vec,
        [&](const Eigen::VectorXd& _d) {
            return barrier->first_derivative(_d[0], dhat);
        },
        fgrad);

    grad << barrier->second_derivative(d, dhat);

    CAPTURE(dhat, d, fgrad(0), grad(0), use_dist_sqr);
    CHECK(fd::compare_gradient(fgrad, grad));
}

TEST_CASE("Physical barrier", "[barrier]")
{
    const bool use_dist_sqr = GENERATE(false, true);

    ipc::ClampedLogBarrier original_barrier;
    PhysicalBarrier new_barrier(use_dist_sqr);

    const double dhat =
        pow(10, GENERATE_COPY(range(use_dist_sqr ? -5 : -5, 0)));

    const double d = GENERATE_COPY(take(10, random(dhat / 2, 0.9 * dhat)));

    const double p_d = (use_dist_sqr ? d : 1) * d;
    const double p_dhat = (use_dist_sqr ? dhat : 1) * dhat;
    const double divisor = original_barrier.units(p_dhat) / dhat;

    CAPTURE(use_dist_sqr, dhat, d, divisor, p_d, p_dhat);

    const double b_original = original_barrier(p_d, p_dhat) / divisor;
    const double b_new = new_barrier(p_d, p_dhat);

    CHECK(b_original == Catch::Approx(b_new));

    const double b_original_first_derivative =
        original_barrier.first_derivative(p_d, p_dhat) / divisor;
    const double b_new_first_derivative =
        new_barrier.first_derivative(p_d, p_dhat);

    CHECK(b_original_first_derivative == Catch::Approx(b_new_first_derivative));

    const double b_original_second_derivative =
        original_barrier.second_derivative(p_d, p_dhat) / divisor;
    const double b_new_second_derivative =
        new_barrier.second_derivative(p_d, p_dhat);

    CHECK(
        b_original_second_derivative == Catch::Approx(b_new_second_derivative));
}