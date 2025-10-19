#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <ipc/distance/point_point.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <ipc/smooth_contact/primitives/edge3.hpp>
#include <ipc/smooth_contact/distance/primitive_distance.hpp>

#include <finitediff.hpp>
#include <igl/PI.h>

using namespace ipc;

namespace {
double edge_edge_distance_stacked(const Eigen::VectorXd& x)
{
    assert(x.size() == 12);
    return edge_edge_distance(
        x.segment<3>(0), x.segment<3>(3), x.segment<3>(6), x.segment<3>(9));
};
} // namespace

TEST_CASE("Edge-edge distance", "[distance][edge-edge]")
{
    double e0y = GENERATE(-10, -1, -1e-4, 0, 1e-4, 1, 10);
    Eigen::Vector3d e00(-1, e0y, 0);
    Eigen::Vector3d e01(1, e0y, 0);
    Eigen::Vector3d e10(0, 0, -1);
    Eigen::Vector3d e11(0, 0, 1);

    Eigen::Vector3d e0_closest, e1_closest;
    double shiftx = GENERATE(-2, 0, 2);
    double shiftz = GENERATE(-2, 0, 2);
    double e0x = shiftx + GENERATE(-1, -0.5, 0, 0.5, 1);
    double e0z = shiftz + GENERATE(-1, -0.5, 0, 0.5, 1);
    e00.x() += e0x;
    e01.x() += e0x;
    e00.z() += e0z;
    e01.z() += e0z;
    e0_closest =
        shiftx > 1 ? e00 : (shiftx < -1 ? e01 : Eigen::Vector3d(0, e0y, e0z));
    e1_closest =
        shiftz > 1 ? e11 : (shiftz < -1 ? e10 : Eigen::Vector3d(0, 0, e0z));

    double distance = edge_edge_distance(e00, e01, e10, e11);
    double expected_distance = point_point_distance(e0_closest, e1_closest);

    CAPTURE(e0x, e0y, e0z, edge_edge_distance_type(e00, e01, e10, e11));
    CHECK(distance == Catch::Approx(expected_distance).margin(1e-12));

    distance =
        PrimitiveDistanceTemplate<Edge3, Edge3, double>::compute_distance(
            (Vector12d() << e00, e01, e10, e11).finished(),
            edge_edge_distance_type(e00, e01, e10, e11));
    CHECK(distance == Catch::Approx(expected_distance));
}

TEST_CASE("Edge-edge distance !EA_EB", "[distance][edge-edge]")
{
    double alpha = GENERATE(range(-1.0, 2.0, 0.1));
    double s = GENERATE(range(-10.0, 10.0, 1.0));
    if (s == 0)
        return;
    const bool swap_ea = GENERATE(false, true);
    const bool swap_eb = GENERATE(false, true);
    const bool swap_edges = GENERATE(false, true);
    const int n_random_edges = 20;

    const auto sign = [](double x) { return x < 0 ? -1 : 1; };

    for (int i = 0; i < n_random_edges; i++) {
        Eigen::Vector3d ea0 = Eigen::Vector3d::Random();
        Eigen::Vector3d ea1 = Eigen::Vector3d::Random();
        Eigen::Vector3d n =
            (ea1 - ea0).cross(Eigen::Vector3d::UnitX()).normalized();

        Eigen::Vector3d eb0 = ((ea1 - ea0) * alpha + ea0) + s * n;
        if (alpha < 0) {
            alpha = 0;
            n = eb0 - ea0;
            s = n.norm();
            n /= s;
        } else if (alpha > 1) {
            alpha = 1;
            n = eb0 - ea1;
            s = n.norm();
            n /= s;
        }
        REQUIRE(s != 0);

        Eigen::Vector3d eb1 =
            ((ea1 - ea0) * alpha + ea0) + (s + sign(s) * 1) * n;

        if (swap_ea) {
            std::swap(ea0, ea1);
        }
        if (swap_eb) {
            std::swap(eb0, eb1);
        }
        if (swap_edges) {
            std::swap(ea0, eb0);
            std::swap(ea1, eb1);
        }

        double distance = edge_edge_distance(ea0, ea1, eb0, eb1);

        CAPTURE(alpha, s, swap_ea, swap_eb, swap_edges);

        CHECK(distance == Catch::Approx(s * s).margin(1e-15));

        distance =
            PrimitiveDistanceTemplate<Edge3, Edge3, double>::compute_distance(
                (Vector12d() << ea0, ea1, eb0, eb1).finished(),
                edge_edge_distance_type(ea0, ea1, eb0, eb1));
        CHECK(distance == Catch::Approx(s * s).margin(1e-15));
    }
}

TEST_CASE("Edge-edge distance EA_EB", "[distance][edge-edge]")
{
    double alpha = GENERATE(take(5, random(0.01, 0.99)));
    double beta = GENERATE(take(5, random(0.01, 0.99)));
    double s = GENERATE(range(-5.0, 5.0, 1.0));
    const bool swap_ea = GENERATE(false, true);
    const bool swap_eb = GENERATE(false, true);
    const bool swap_edges = GENERATE(false, true);
    const int n_random_edges = 20;

    for (int i = 0; i < n_random_edges; i++) {
        Eigen::Vector3d ea1 = Eigen::Vector3d::Random();
        Eigen::Vector3d eb1 = Eigen::Vector3d::Random();

        Eigen::Vector3d ea0 = alpha / (alpha - 1) * ea1;
        Eigen::Vector3d eb0 = beta / (beta - 1) * eb1;

        Eigen::Vector3d n = ea1.cross(eb1);
        REQUIRE(n.norm() != 0);
        n.normalize();

        ea0 += -s / 2 * n;
        ea1 += -s / 2 * n;
        eb0 += s / 2 * n;
        eb1 += s / 2 * n;

        if (swap_ea) {
            std::swap(ea0, ea1);
        }
        if (swap_eb) {
            std::swap(eb0, eb1);
        }
        if (swap_edges) {
            std::swap(ea0, eb0);
            std::swap(ea1, eb1);
        }

        double distance = edge_edge_distance(ea0, ea1, eb0, eb1);
        CHECK(distance == Catch::Approx(s * s).margin(1e-15));

        distance =
            PrimitiveDistanceTemplate<Edge3, Edge3, double>::compute_distance(
                (Vector12d() << ea0, ea1, eb0, eb1).finished(),
                edge_edge_distance_type(ea0, ea1, eb0, eb1));
        CHECK(distance == Catch::Approx(s * s).margin(1e-15));
    }
}

TEST_CASE("Edge-edge distance parallel", "[distance][edge-edge][parallel]")
{
    double alpha = GENERATE(take(10, random(0.01, 0.99)));
    double s = GENERATE(take(10, random(-5.0, 5.0)));
    const int n_random_edges = 20;

    for (int i = 0; i < n_random_edges; i++) {
        const Eigen::Vector3d ea0 = Eigen::Vector3d::Random();
        const Eigen::Vector3d ea1 = Eigen::Vector3d::Random();
        const double edge_len = (ea1 - ea0).norm();

        Eigen::Vector3d n = (ea1 - ea0).cross(Eigen::Vector3d::UnitX());
        REQUIRE(n.norm() != 0);
        n.normalize();

        const Eigen::Vector3d eb0 = (ea1 - ea0) * alpha + ea0 + s * n;
        const Eigen::Vector3d eb1 = (ea1 - ea0) + eb0;
        REQUIRE(
            (ea1 - ea0).dot(eb1 - eb0) == Catch::Approx(edge_len * edge_len));
        REQUIRE(
            (ea1 - ea0).cross(eb1 - eb0).norm()
            == Catch::Approx(0).margin(1e-14));

        const double distance = edge_edge_distance(ea0, ea1, eb0, eb1);
        CHECK(distance == Catch::Approx(s * s).margin(1e-15));

        for (int dtype = 0; dtype < int(EdgeEdgeDistanceType::EA_EB); dtype++) {
            const double distance2 = edge_edge_distance(
                ea0, ea1, eb0, eb1, EdgeEdgeDistanceType(dtype));
            CAPTURE(dtype);
            CHECK(distance <= Catch::Approx(distance2));
        }
    }
}

TEST_CASE("Edge-edge distance degenerate case", "[distance][edge-edge]")
{
    double e0y = GENERATE(-10, -1, -1e-4, 0, 1e-4, 1, 10);
    Eigen::Vector3d e00(-1, e0y, 0);
    Eigen::Vector3d e01(1, e0y, 0);
    Eigen::Vector3d e10(0, 0, -1);
    Eigen::Vector3d e11(0, 0, 1);

    double theta =
        GENERATE(-2, -1.5, -1, -0.123124, 0, 0.2342352, 0.5, 1, 1.5, 2, 50, 51)
        * igl::PI;
    Eigen::Matrix3d R =
        Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitY()).toRotationMatrix();

    SECTION("e0 rotating")
    {
        e00 = R * e00;
        e01 = R * e01;
    }
    SECTION("e1 rotating")
    {
        e10 = R * e10;
        e11 = R * e11;
    }

    double distance = edge_edge_distance(e00, e01, e10, e11);
    CHECK(distance == Catch::Approx(e0y * e0y).margin(1e-12));
}

TEST_CASE(
    "Edge-edge distance degenerate case not overlapping",
    "[distance][edge-edge]")
{
    double e0y = GENERATE(-10, -1, -1e-4, 0, 1e-4, 1, 10);
    double gap = GENERATE(0, 0.01, 0.1, 1);
    Eigen::Vector3d e00(gap, e0y, 0);
    Eigen::Vector3d e01(1, e0y, 0);
    Eigen::Vector3d e10(-1, 0, 0);
    Eigen::Vector3d e11(-gap, 0, 0);

    SECTION("original order") { }
    SECTION("swap e0") { std::swap(e00, e01); }
    SECTION("swap e1") { std::swap(e10, e11); }
    SECTION("swap e0 and e1")
    {
        std::swap(e00, e01);
        std::swap(e10, e11);
    }

    double distance = edge_edge_distance(e00, e01, e10, e11);
    double expected_distance = point_point_distance(
        Eigen::Vector3d(gap, e0y, 0), Eigen::Vector3d(-gap, 0, 0));

    CHECK(distance == Catch::Approx(expected_distance).margin(1e-12));
}

TEST_CASE("Edge-edge distance gradient", "[distance][edge-edge][gradient]")
{
    double e0y = GENERATE(-10, -1, -1e-4, 0, 1e-4, 1, 10);
    Eigen::Vector3d e00(-1, e0y, 0);
    Eigen::Vector3d e01(1, e0y, 0);
    Eigen::Vector3d e10(0, 0, -1);
    Eigen::Vector3d e11(0, 0, 1);

    Eigen::Vector3d e0_closest, e1_closest;
    double shiftx = GENERATE(-2, 0, 2);
    double shiftz = GENERATE(-2, 0, 2);
    double e0x = shiftx + GENERATE(-1, -0.5, 0, 0.5, 1);
    double e0z = shiftz + GENERATE(-1, -0.5, 0, 0.5, 1);
    e00.x() += e0x;
    e01.x() += e0x;
    e00.z() += e0z;
    e01.z() += e0z;
    e0_closest =
        shiftx > 1 ? e00 : (shiftx < -1 ? e01 : Eigen::Vector3d(0, e0y, e0z));
    e1_closest =
        shiftz > 1 ? e11 : (shiftz < -1 ? e10 : Eigen::Vector3d(0, 0, e0z));

    const Vector12d grad = edge_edge_distance_gradient(e00, e01, e10, e11);

    Vector12d x;
    x << e00, e01, e10, e11;

    // Compute the gradient using finite differences
    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, edge_edge_distance_stacked, fgrad);

    CAPTURE(e0x, e0y, e0z, edge_edge_distance_type(e00, e01, e10, e11));
    CHECK(fd::compare_gradient(grad, fgrad));

    Vector<ADGrad<12>, 12> X1 = slice_positions<ADGrad<12>, 12, 1>(
        (Vector12d() << e00, e01, e10, e11).finished());
    ADGrad<12> dist1 =
        PrimitiveDistanceTemplate<Edge3, Edge3, ADGrad<12>>::compute_distance(
            X1, edge_edge_distance_type(e00, e01, e10, e11));
    CHECK(fd::compare_gradient(dist1.grad, fgrad));
}

TEST_CASE(
    "Parallel edge-edge distance gradient", "[distance][edge-edge][gradient]")
{
    // Generate a geometric space of
    double angle = 0;
    SECTION("Almost parallel")
    {
        double exponent = GENERATE(range(-6, 3));
        angle = pow(10, exponent) * igl::PI / 180.0;
    }
    // SECTION("Parallel") { angle = 0; }

    Eigen::Vector3d e00(-1.0, 0, 0), e01(1.0, 0, 0),
        e10(cos(angle), 1, sin(angle)),
        e11(cos(angle + igl::PI), 1, sin(angle + igl::PI));

    double distance = edge_edge_distance(e00, e01, e10, e11);
    const Vector12d grad = edge_edge_distance_gradient(e00, e01, e10, e11);

    Vector12d x;
    x << e00, e01, e10, e11;

    // Compute the gradient using finite differences
    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, edge_edge_distance_stacked, fgrad);

    CAPTURE(angle, (grad - fgrad).squaredNorm());
    CHECK(distance == Catch::Approx(1.0));
    CHECK(fd::compare_gradient(grad, fgrad));
    // CHECK(distance.Hess.squaredNorm() != Catch::Approx(0.0));
}

TEST_CASE(
    "Edge-edge closest direction gradient !EA_EB", "[distance][edge-edge]")
{
    double alpha = GENERATE(range(-1.0, 2.0, 0.1));
    double s = GENERATE(range(-10.0, 10.0, 1.0));
    if (s == 0)
        return;
    const bool swap_ea = GENERATE(false, true);
    const bool swap_eb = GENERATE(false, true);
    const bool swap_edges = GENERATE(false, true);
    const int n_random_edges = 20;

    const auto sign = [](double x) { return x < 0 ? -1 : 1; };

    for (int i = 0; i < n_random_edges; i++) {
        Eigen::Vector3d ea0 = Eigen::Vector3d::Random();
        Eigen::Vector3d ea1 = Eigen::Vector3d::Random();
        Eigen::Vector3d n =
            (ea1 - ea0).cross(Eigen::Vector3d::UnitX()).normalized();

        Eigen::Vector3d eb0 = ((ea1 - ea0) * alpha + ea0) + s * n;
        if (alpha < 0) {
            alpha = 0;
            n = eb0 - ea0;
            s = n.norm();
            n /= s;
        } else if (alpha > 1) {
            alpha = 1;
            n = eb0 - ea1;
            s = n.norm();
            n /= s;
        }
        REQUIRE(s != 0);

        Eigen::Vector3d eb1 =
            ((ea1 - ea0) * alpha + ea0) + (s + sign(s) * 1) * n;

        if (swap_ea) {
            std::swap(ea0, ea1);
        }
        if (swap_eb) {
            std::swap(eb0, eb1);
        }
        if (swap_edges) {
            std::swap(ea0, eb0);
            std::swap(ea1, eb1);
        }

        double distance = edge_edge_distance(ea0, ea1, eb0, eb1);

        CAPTURE(alpha, s, swap_ea, swap_eb, swap_edges);

        CHECK(distance == Catch::Approx(s * s).margin(1e-15));

        Vector12d x;
        x << ea0, ea1, eb0, eb1;
        const auto dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);

        distance =
            PrimitiveDistanceTemplate<Edge3, Edge3, double>::compute_distance(
                x, dtype);
        CHECK(distance == Catch::Approx(s * s).margin(1e-15));

        Vector<ADGrad<12>, 12> X = slice_positions<ADGrad<12>, 12, 1>(x);
        Vector<ADGrad<12>, 3> direc = PrimitiveDistanceTemplate<
            Edge3, Edge3, ADGrad<12>>::compute_closest_direction(X, dtype);
        CHECK(direc.squaredNorm().val == Catch::Approx(s * s).margin(1e-15));

        // Compute the gradient using finite differences
        Eigen::MatrixXd fjac;
        fd::finite_jacobian(
            x,
            [dtype](const Eigen::VectorXd& _x) {
                return PrimitiveDistanceTemplate<
                    Edge3, Edge3, double>::compute_closest_direction(_x, dtype);
            },
            fjac);

        for (int d = 0; d < 3; d++) {
            CHECK(
                fd::compare_jacobian(
                    fjac.row(d), direc(d).grad.transpose(), 1e-5));
        }
    }
}

TEST_CASE("Edge normal term", "[distance][edge-edge][gradient]")
{
    OrientationTypes otypes;
    otypes.set_size(2);
    const double alpha = 0.85;
    const double beta = 0.2;
    SmoothContactParameters params { 1e-3, 1, 0, alpha, beta, 2 };

    ipc::Vector3d dn, e0, e1, f0, f1;
    e0 << 0, 0, 0;
    e1 << 1, 0, 0;
    dn << 0, 0, 1;
    f0 << 0.4, 0.3, GENERATE(take(10, random(-0.2, 0.2)));
    f1 << 0.6, -0.2, GENERATE(take(10, random(-0.2, 0.2)));

    ipc::Vector15d x;
    x << dn, e0, e1, f0, f1;

    const auto [val, grad, hess] = smooth_edge3_normal_term_hessian(
        dn, e0, e1, f0, f1, alpha, beta, otypes);

    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        x,
        [&](const Eigen::VectorXd& y) {
            return smooth_edge3_normal_term(
                y.head(3), y.segment(3, 3), y.segment(6, 3), y.segment(9, 3),
                y.segment(12, 3), alpha, beta, otypes);
        },
        fgrad, fd::AccuracyOrder::FOURTH, 1e-5);

    CHECK((fgrad - grad).norm() / 1e-8 <= grad.norm());

    Eigen::MatrixXd fhess;
    fd::finite_jacobian(
        x,
        [&](const Eigen::VectorXd& y) {
            return std::get<1>(smooth_edge3_normal_term_gradient(
                y.head(3), y.segment(3, 3), y.segment(6, 3), y.segment(9, 3),
                y.segment(12, 3), alpha, beta, otypes));
        },
        fhess, fd::AccuracyOrder::SECOND, 1e-7);

    CHECK((fhess - hess).norm() / 1e-6 <= hess.norm());
}

TEMPLATE_TEST_CASE_SIG(
    "Edge-edge templated functions !EA_EB",
    "[distance][edge-edge]",
    ((int ndofs), ndofs),
    9,
    12,
    13)
{
    double alpha = GENERATE(range(-1.0, 2.0, 0.4));
    double s = GENERATE(range(-10.0, 10.0, 5.0));
    if (s == 0)
        return;
    const bool swap_ea = GENERATE(false, true);
    const bool swap_eb = GENERATE(false, true);
    const bool swap_edges = GENERATE(false, true);
    const int n_random_edges = 5;

    const auto sign = [](double x) { return x < 0 ? -1 : 1; };

    for (int i = 0; i < n_random_edges; i++) {
        Eigen::Vector3d ea0 = Eigen::Vector3d::Random();
        Eigen::Vector3d ea1 = Eigen::Vector3d::Random();
        Eigen::Vector3d n =
            (ea1 - ea0).cross(Eigen::Vector3d::UnitX()).normalized();

        Eigen::Vector3d eb0 = ((ea1 - ea0) * alpha + ea0) + s * n;
        if (alpha < 0) {
            alpha = 0;
            n = eb0 - ea0;
            s = n.norm();
            n /= s;
        } else if (alpha > 1) {
            alpha = 1;
            n = eb0 - ea1;
            s = n.norm();
            n /= s;
        }
        REQUIRE(s != 0);

        Eigen::Vector3d eb1 =
            ((ea1 - ea0) * alpha + ea0) + (s + sign(s) * 1) * n;

        if (swap_ea) {
            std::swap(ea0, ea1);
        }
        if (swap_eb) {
            std::swap(eb0, eb1);
        }
        if (swap_edges) {
            std::swap(ea0, eb0);
            std::swap(ea1, eb1);
        }

        Vector12d x;
        x << ea0, ea1, eb0, eb1;
        const auto dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);

        auto distance1 = PrimitiveDistanceTemplate<
            Edge3, Edge3, ADGrad<ndofs>>::compute_distance(x, dtype);
        CHECK(distance1.val == Catch::Approx(s * s).margin(1e-15));

        auto distance2 = PrimitiveDistanceTemplate<
            Edge3, Edge3, ADHessian<ndofs>>::compute_distance(x, dtype);
        CHECK(distance2.val == Catch::Approx(s * s).margin(1e-15));
    }
}

TEMPLATE_TEST_CASE_SIG(
    "Edge-edge templated functions EA_EB",
    "[distance][edge-edge]",
    ((int ndofs), ndofs),
    9,
    12,
    13)
{
    Eigen::Vector3d e00(-1, 0, 0);
    Eigen::Vector3d e01(1, 0, 0);
    Eigen::Vector3d e10(0, 0, -1);
    Eigen::Vector3d e11(0, 0, 1);

    Vector12d x;
    x << e00, e01, e10, e11;
    const auto dtype = edge_edge_distance_type(e00, e01, e10, e11);
    double expected_distance = edge_edge_distance(e00, e01, e10, e11, dtype);

    auto distance1 = PrimitiveDistanceTemplate<
        Edge3, Edge3, ADGrad<ndofs>>::compute_distance(x, dtype);
    CHECK(distance1.val == Catch::Approx(expected_distance).margin(1e-15));

    auto distance2 = PrimitiveDistanceTemplate<
        Edge3, Edge3, ADHessian<ndofs>>::compute_distance(x, dtype);
    CHECK(distance2.val == Catch::Approx(expected_distance).margin(1e-15));

    auto direc1 = PrimitiveDistanceTemplate<
        Edge3, Edge3, ADGrad<ndofs>>::compute_closest_direction(x, dtype);
    CHECK(
        direc1.squaredNorm().val
        == Catch::Approx(expected_distance).margin(1e-15));

    auto direc2 = PrimitiveDistanceTemplate<
        Edge3, Edge3, ADHessian<ndofs>>::compute_closest_direction(x, dtype);
    CHECK(
        direc2.squaredNorm().val
        == Catch::Approx(expected_distance).margin(1e-15));
}