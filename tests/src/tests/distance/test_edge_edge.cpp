#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <ipc/collision_mesh.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/smooth_contact/primitives/edge3.hpp>
#include <ipc/smooth_contact/distance/primitive_distance.hpp>
#include <ipc/smooth_contact/distance/point_edge.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <finitediff.hpp>
#include <igl/PI.h>
#include <igl/edges.h>

using namespace ipc;

namespace {
double edge_edge_distance_stacked(const Eigen::VectorXd& x)
{
    assert(x.size() == 12);
    return edge_edge_distance(
        x.segment<3>(0), x.segment<3>(3), x.segment<3>(6), x.segment<3>(9));
}
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

    Eigen::Vector<ADGrad<12>, 12> X1 = slice_positions<ADGrad<12>, 12, 1>(
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

        Eigen::Vector<ADGrad<12>, 12> X = slice_positions<ADGrad<12>, 12, 1>(x);
        Eigen::Vector<ADGrad<12>, 3> direc = PrimitiveDistanceTemplate<
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
    const double alpha = 0.85;
    const double beta = 0.2;
    SmoothContactParameters params { 1e-3, 1, 0, alpha, beta, 2 };

    Eigen::Vector3d dn, e0, e1, f0, f1;
    e0 << 0, 0, 0;
    e1 << 1, 0, 0;
    dn << 0, 0, 1;
    f0 << 0.4, 0.3, GENERATE(take(10, random(-0.2, 0.2)));
    f1 << 0.6, -0.2, GENERATE(take(10, random(-0.2, 0.2)));

    // Build a minimal mesh: edge (0,1) with two adjacent faces sharing it.
    // Vertices: 0=e0, 1=e1, 2=f0, 3=f1
    Eigen::MatrixXd V(4, 3);
    V.row(0) = e0.transpose();
    V.row(1) = e1.transpose();
    V.row(2) = f0.transpose();
    V.row(3) = f1.transpose();

    // Two faces sharing edge (0,1), wound so the edge is traversed in
    // opposite directions (proper orientation).
    Eigen::MatrixXi F(2, 3);
    F << 2, 0, 1, 3, 1, 0;

    Eigen::MatrixXi E;
    igl::edges(F, E);

    // Find the edge index for the edge (0,1)
    int edge_id = -1;
    for (int i = 0; i < E.rows(); i++) {
        if ((E(i, 0) == 0 && E(i, 1) == 1) || (E(i, 0) == 1 && E(i, 1) == 0)) {
            edge_id = i;
            break;
        }
    }
    REQUIRE(edge_id >= 0);

    // Construct CollisionMesh with all vertices marked as orientable
    std::vector<bool> include_vertex(4, true);
    std::vector<bool> orient_vertex(4, true);
    CollisionMesh mesh(include_vertex, orient_vertex, V, E, F);

    // Construct Edge3 from the mesh (this sets otypes, faces, etc.)
    auto edge3 = std::make_unique<Edge3>(edge_id, mesh, V, dn, params);

    // X layout for the member functions: rows are [e0, e1, f0, f1]
    // matching m_vertex_ids order from the Edge3 constructor.
    Eigen::MatrixX3d X(4, 3);
    // The Edge3 constructor orders face-opposite vertices by adjacency order.
    // We reconstruct X from the Edge3's vertex_ids to match internal ordering.
    const auto& vids = edge3->vertex_ids();
    for (int i = 0; i < 4; i++) {
        X.row(i) = V.row(vids[i]);
    }

    // Flatten [dn, X] into a single vector for finite differences
    Eigen::VectorXd x(3 + 4 * 3);
    x.head<3>() = dn;
    for (int i = 0; i < 4; i++) {
        x.segment<3>(3 + i * 3) = X.row(i);
    }

    const auto [val, grad, hess] =
        edge3->smooth_edge3_normal_term_hessian(dn, X, alpha, beta);

    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        x,
        [&](const Eigen::VectorXd& y) {
            Eigen::MatrixX3d Y = fd::unflatten(y.tail(12), 3);
            return std::get<0>(edge3->smooth_edge3_normal_term_gradient(
                y.head(3), Y, alpha, beta));
        },
        fgrad, fd::AccuracyOrder::FOURTH, 1e-5);

    CHECK((fgrad - grad).norm() / 1e-8 <= grad.norm());

    Eigen::MatrixXd fhess;
    fd::finite_jacobian(
        x,
        [&](const Eigen::VectorXd& y) -> Eigen::VectorXd {
            Eigen::MatrixX3d Y = fd::unflatten(y.tail(12), 3);
            return std::get<1>(edge3->smooth_edge3_normal_term_gradient(
                y.head(3), Y, alpha, beta));
        },
        fhess, fd::AccuracyOrder::SECOND, 1e-7);

    CHECK((fhess - hess).norm() / 1e-6 <= hess.norm());
}

// =============================================================================
// Helper: build a mesh with an edge having the given number of face neighbors,
// construct an Edge3, and return the edge3, mesh, edge_id, and vertex-id map.
// =============================================================================
namespace {
struct Edge3TestFixture {
    std::unique_ptr<Edge3> edge3;
    CollisionMesh mesh;
    int edge_id;

    /// Build a fixture for an edge along the x-axis with nf face neighbors.
    /// orient: if true, vertices are marked orientable (normal term active).
    /// Vertices: 0=e0=(0,0,0), 1=e1=(1,0,0), 2..2+nf-1 = opposite vertices.
    /// Faces wound so the shared edge is traversed in opposite directions in
    /// adjacent faces: face a = (2+a, a%2==0 ? 0 : 1, a%2==0 ? 1 : 0).
    static Edge3TestFixture build(
        const Eigen::Vector3d& e0,
        const Eigen::Vector3d& e1,
        const Eigen::MatrixX3d& face_verts, // nf x 3
        const Eigen::Vector3d& dn,
        const SmoothContactParameters& params,
        bool orient)
    {
        const int nf = static_cast<int>(face_verts.rows());
        const int nv = 2 + nf;

        Eigen::MatrixXd V(nv, 3);
        V.row(0) = e0.transpose();
        V.row(1) = e1.transpose();
        for (int a = 0; a < nf; a++) {
            V.row(2 + a) = face_verts.row(a);
        }

        // Wind faces so edge 0-1 appears in alternating directions.
        Eigen::MatrixXi F(nf, 3);
        for (int a = 0; a < nf; a++) {
            if (a % 2 == 0) {
                F.row(a) << 2 + a, 0, 1;
            } else {
                F.row(a) << 2 + a, 1, 0;
            }
        }

        Eigen::MatrixXi E;
        igl::edges(F, E);

        int eid = -1;
        for (int i = 0; i < E.rows(); i++) {
            if ((E(i, 0) == 0 && E(i, 1) == 1)
                || (E(i, 0) == 1 && E(i, 1) == 0)) {
                eid = i;
                break;
            }
        }

        std::vector<bool> include_vertex(nv, true);
        std::vector<bool> orient_vertex(nv, orient);
        CollisionMesh m(include_vertex, orient_vertex, V, E, F);

        auto e3 = std::make_unique<Edge3>(eid, m, V, dn, params);

        return { std::move(e3), std::move(m), eid };
    }

    /// Assemble X matrix from the Edge3's vertex_ids and given vertex
    /// positions.
    Eigen::MatrixX3d build_X(const Eigen::MatrixXd& V) const
    {
        const auto& vids = edge3->vertex_ids();
        Eigen::MatrixX3d X(static_cast<int>(vids.size()), 3);
        for (int i = 0; i < static_cast<int>(vids.size()); i++) {
            X.row(i) = V.row(vids[i]);
        }
        return X;
    }

    /// Assemble tangent directions from X (for tangent term tests).
    Eigen::MatrixX3d build_tangents(const Eigen::MatrixX3d& X) const
    {
        const int nn = edge3->n_face_neighbors();
        Eigen::MatrixX3d tangents(nn, 3);
        for (int a = 0; a < nn; a++) {
            tangents.row(a) = PointEdgeDistance<double, 3>::
                point_line_closest_point_direction(
                    X.row(2 + a), X.row(0), X.row(1));
        }
        return tangents;
    }

    /// Flatten [dn, tangents] or [dn, X] into a single FD vector.
    static Eigen::VectorXd
    flatten_with_dir(const Eigen::Vector3d& dn, const Eigen::MatrixX3d& M)
    {
        const int rows = static_cast<int>(M.rows());
        Eigen::VectorXd x(3 + rows * 3);
        x.head<3>() = dn;
        for (int i = 0; i < rows; i++) {
            x.segment<3>(3 + i * 3) = M.row(i);
        }
        return x;
    }
};
} // namespace

// =============================================================================
// Tangent term hessian — 2 neighbors, both VARIANT
// =============================================================================
TEST_CASE(
    "Edge tangent term hessian 2 neighbors", "[distance][edge-edge][gradient]")
{
    // alpha_t=1, beta_t=0 ⇒ VARIANT when tangent check val ∈ (-1, 0),
    // i.e. dn·t / |t| ∈ (0, 1).
    // With dn=(0,0,1), t=(0, y, z), we need z/|t| ∈ (0, 1) ⇒ z > 0.
    // alpha_n, beta_n chosen so normal term stays ONE (no interaction).
    const double alpha_t = 1.0, beta_t = 0.0;
    SmoothContactParameters params(1e-3, alpha_t, beta_t, 0.85, 0.2, 2);

    Eigen::Vector3d e0(0, 0, 0), e1(1, 0, 0), dn(0, 0, 1);
    // Both f0, f1 above the edge so tangent z > 0 ⇒ VARIANT.
    Eigen::Vector3d f0(0.4, 0.3, GENERATE(take(5, random(0.02, 0.15))));
    Eigen::Vector3d f1(0.6, -0.2, GENERATE(take(5, random(0.02, 0.15))));

    Eigen::MatrixX3d fv(2, 3);
    fv.row(0) = f0.transpose();
    fv.row(1) = f1.transpose();

    auto fix = Edge3TestFixture::build(e0, e1, fv, dn, params, true);
    REQUIRE(fix.edge3->n_face_neighbors() == 2);

    Eigen::MatrixXd V(4, 3);
    V.row(0) = e0.transpose();
    V.row(1) = e1.transpose();
    V.row(2) = f0.transpose();
    V.row(3) = f1.transpose();

    Eigen::MatrixX3d X = fix.build_X(V);
    Eigen::MatrixX3d tangents = fix.build_tangents(X);

    // Evaluate hessian
    const auto [val, grad, hess] = fix.edge3->smooth_edge3_tangent_term_hessian(
        dn, tangents, alpha_t, beta_t);

    // FD vector: [dn(3), tangent_0(3), tangent_1(3)]
    Eigen::VectorXd x = Edge3TestFixture::flatten_with_dir(dn, tangents);

    // Gradient via FD
    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        x,
        [&](const Eigen::VectorXd& y) {
            Eigen::MatrixX3d T = fd::unflatten(y.tail(y.size() - 3), 3);
            return std::get<0>(fix.edge3->smooth_edge3_tangent_term_gradient(
                y.head(3), T, alpha_t, beta_t));
        },
        fgrad, fd::AccuracyOrder::FOURTH, 1e-5);

    CHECK((fgrad - grad).norm() / 1e-8 <= grad.norm());

    // Hessian via FD of gradient
    Eigen::MatrixXd fhess;
    fd::finite_jacobian(
        x,
        [&](const Eigen::VectorXd& y) -> Eigen::VectorXd {
            Eigen::MatrixX3d T = fd::unflatten(y.tail(y.size() - 3), 3);
            return std::get<1>(fix.edge3->smooth_edge3_tangent_term_gradient(
                y.head(3), T, alpha_t, beta_t));
        },
        fhess, fd::AccuracyOrder::SECOND, 1e-7);

    CHECK((fhess - hess).norm() / 1e-6 <= hess.norm());
}

// =============================================================================
// Tangent term hessian — 1 neighbor (boundary edge)
// =============================================================================
TEST_CASE(
    "Edge tangent term hessian 1 neighbor", "[distance][edge-edge][gradient]")
{
    const double alpha_t = 1.0, beta_t = 0.0;
    SmoothContactParameters params(1e-3, alpha_t, beta_t, 0.85, 0.2, 2);

    Eigen::Vector3d e0(0, 0, 0), e1(1, 0, 0), dn(0, 0, 1);
    Eigen::Vector3d f0(0.5, 0.4, GENERATE(take(5, random(0.02, 0.15))));

    Eigen::MatrixX3d fv(1, 3);
    fv.row(0) = f0.transpose();

    auto fix = Edge3TestFixture::build(e0, e1, fv, dn, params, true);
    REQUIRE(fix.edge3->n_face_neighbors() == 1);

    Eigen::MatrixXd V(3, 3);
    V.row(0) = e0.transpose();
    V.row(1) = e1.transpose();
    V.row(2) = f0.transpose();

    Eigen::MatrixX3d X = fix.build_X(V);
    Eigen::MatrixX3d tangents = fix.build_tangents(X);

    const auto [val, grad, hess] = fix.edge3->smooth_edge3_tangent_term_hessian(
        dn, tangents, alpha_t, beta_t);

    Eigen::VectorXd x = Edge3TestFixture::flatten_with_dir(dn, tangents);

    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        x,
        [&](const Eigen::VectorXd& y) {
            Eigen::MatrixX3d T = fd::unflatten(y.tail(y.size() - 3), 3);
            return std::get<0>(fix.edge3->smooth_edge3_tangent_term_gradient(
                y.head(3), T, alpha_t, beta_t));
        },
        fgrad, fd::AccuracyOrder::FOURTH, 1e-5);

    CHECK((fgrad - grad).norm() / 1e-8 <= grad.norm());

    Eigen::MatrixXd fhess;
    fd::finite_jacobian(
        x,
        [&](const Eigen::VectorXd& y) -> Eigen::VectorXd {
            Eigen::MatrixX3d T = fd::unflatten(y.tail(y.size() - 3), 3);
            return std::get<1>(fix.edge3->smooth_edge3_tangent_term_gradient(
                y.head(3), T, alpha_t, beta_t));
        },
        fhess, fd::AccuracyOrder::SECOND, 1e-7);

    CHECK((fhess - hess).norm() / 1e-6 <= hess.norm());
}

// =============================================================================
// Normal term hessian — non-orientable edge (early-return path)
// =============================================================================
TEST_CASE("Edge normal term non-orientable", "[distance][edge-edge][gradient]")
{
    const double alpha_n = 0.85, beta_n = 0.2;
    SmoothContactParameters params(1e-3, 1, 0, alpha_n, beta_n, 2);

    Eigen::Vector3d e0(0, 0, 0), e1(1, 0, 0), dn(0, 0, 1);
    Eigen::Vector3d f0(0.4, 0.3, 0.1), f1(0.6, -0.2, 0.1);

    Eigen::MatrixX3d fv(2, 3);
    fv.row(0) = f0.transpose();
    fv.row(1) = f1.transpose();

    // orient=false ⇒ orientable will be false inside Edge3
    auto fix = Edge3TestFixture::build(e0, e1, fv, dn, params, false);
    REQUIRE(fix.edge3->n_face_neighbors() == 2);

    Eigen::MatrixXd V(4, 3);
    V.row(0) = e0.transpose();
    V.row(1) = e1.transpose();
    V.row(2) = f0.transpose();
    V.row(3) = f1.transpose();

    Eigen::MatrixX3d X = fix.build_X(V);

    const auto [val, grad, hess] =
        fix.edge3->smooth_edge3_normal_term_hessian(dn, X, alpha_n, beta_n);

    // Non-orientable ⇒ trivially 1 with zero derivatives.
    CHECK(val == Catch::Approx(1.0));
    CHECK(grad.norm() == Catch::Approx(0.0).margin(1e-15));
    CHECK(hess.norm() == Catch::Approx(0.0).margin(1e-15));
}

// =============================================================================
// Normal term hessian — all normal types are ONE (early-return path)
// =============================================================================
TEST_CASE("Edge normal term ONE type", "[distance][edge-edge][gradient]")
{
    // Choose alpha_n, beta_n so that dn clearly aligns with face normals
    // and normal_sum >= 1 at construction, forcing all normal types to ONE.
    // With e0=(0,0,0), e1=(1,0,0), dn=(0,0,1), and faces that have normals
    // pointing in the +z direction, the dot product will be ~1.
    const double alpha_n = 0.5, beta_n = 0.1;
    SmoothContactParameters params(1e-3, 1, 0, alpha_n, beta_n, 2);

    Eigen::Vector3d e0(0, 0, 0), e1(1, 0, 0), dn(0, 0, 1);
    // Place face vertices far from z=0 in the +y direction, so that the
    // face normals point nearly purely in +z (aligned with dn).
    Eigen::Vector3d f0(0.4, 0.5, 0.0), f1(0.6, -0.5, 0.0);

    Eigen::MatrixX3d fv(2, 3);
    fv.row(0) = f0.transpose();
    fv.row(1) = f1.transpose();

    auto fix = Edge3TestFixture::build(e0, e1, fv, dn, params, true);
    REQUIRE(fix.edge3->n_face_neighbors() == 2);

    Eigen::MatrixXd V(4, 3);
    V.row(0) = e0.transpose();
    V.row(1) = e1.transpose();
    V.row(2) = f0.transpose();
    V.row(3) = f1.transpose();

    Eigen::MatrixX3d X = fix.build_X(V);

    const auto [val, grad, hess] =
        fix.edge3->smooth_edge3_normal_term_hessian(dn, X, alpha_n, beta_n);

    // normal_type(0) == ONE ⇒ trivially 1 with zero derivatives.
    CHECK(val == Catch::Approx(1.0));
    CHECK(grad.norm() == Catch::Approx(0.0).margin(1e-15));
    CHECK(hess.norm() == Catch::Approx(0.0).margin(1e-15));
}

// =============================================================================
// Normal term hessian — single face neighbor, VARIANT
// =============================================================================
TEST_CASE("Edge normal term single neighbor", "[distance][edge-edge][gradient]")
{
    const double alpha_n = 0.85, beta_n = 0.2;
    SmoothContactParameters params(1e-3, 1, 0, alpha_n, beta_n, 2);

    Eigen::Vector3d e0(0, 0, 0), e1(1, 0, 0), dn(0, 0, 1);
    Eigen::Vector3d f0(0.4, 0.3, GENERATE(take(5, random(-0.2, 0.2))));

    Eigen::MatrixX3d fv(1, 3);
    fv.row(0) = f0.transpose();

    auto fix = Edge3TestFixture::build(e0, e1, fv, dn, params, true);
    REQUIRE(fix.edge3->n_face_neighbors() == 1);

    Eigen::MatrixXd V(3, 3);
    V.row(0) = e0.transpose();
    V.row(1) = e1.transpose();
    V.row(2) = f0.transpose();

    Eigen::MatrixX3d X = fix.build_X(V);

    Eigen::VectorXd x(3 + 3 * 3);
    x.head<3>() = dn;
    for (int i = 0; i < 3; i++) {
        x.segment<3>(3 + i * 3) = X.row(i);
    }

    const auto [val, grad, hess] =
        fix.edge3->smooth_edge3_normal_term_hessian(dn, X, alpha_n, beta_n);

    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        x,
        [&](const Eigen::VectorXd& y) {
            Eigen::MatrixX3d Y = fd::unflatten(y.tail(9), 3);
            return std::get<0>(fix.edge3->smooth_edge3_normal_term_gradient(
                y.head(3), Y, alpha_n, beta_n));
        },
        fgrad, fd::AccuracyOrder::FOURTH, 1e-5);

    CHECK((fgrad - grad).norm() / 1e-8 <= grad.norm());

    Eigen::MatrixXd fhess;
    fd::finite_jacobian(
        x,
        [&](const Eigen::VectorXd& y) -> Eigen::VectorXd {
            Eigen::MatrixX3d Y = fd::unflatten(y.tail(9), 3);
            return std::get<1>(fix.edge3->smooth_edge3_normal_term_gradient(
                y.head(3), Y, alpha_n, beta_n));
        },
        fhess, fd::AccuracyOrder::SECOND, 1e-7);

    CHECK((fhess - hess).norm() / 1e-6 <= hess.norm());
}

// =============================================================================
// Tangent term hessian — mixed types (ONE + VARIANT among neighbors)
// =============================================================================
TEST_CASE(
    "Edge tangent term hessian mixed types", "[distance][edge-edge][gradient]")
{
    // alpha_t=1, beta_t=0.
    // VARIANT when -dn·t/|t| ∈ (-1, 0), i.e. dn·t/|t| ∈ (0, 1) ⇒ z > 0.
    // ONE     when -dn·t/|t| >= 0,       i.e. dn·t/|t| <= 0   ⇒ z <= 0.
    // Place f0 with z > 0 (VARIANT) and f1 with z < 0 (ONE).
    const double alpha_t = 1.0, beta_t = 0.0;
    SmoothContactParameters params(1e-3, alpha_t, beta_t, 0.85, 0.2, 2);

    Eigen::Vector3d e0(0, 0, 0), e1(1, 0, 0), dn(0, 0, 1);
    Eigen::Vector3d f0(0.4, 0.3, GENERATE(take(5, random(0.02, 0.15))));
    Eigen::Vector3d f1(0.6, -0.3, GENERATE(take(5, random(-0.5, -0.1))));

    Eigen::MatrixX3d fv(2, 3);
    fv.row(0) = f0.transpose();
    fv.row(1) = f1.transpose();

    auto fix = Edge3TestFixture::build(e0, e1, fv, dn, params, true);
    REQUIRE(fix.edge3->n_face_neighbors() == 2);

    Eigen::MatrixXd V(4, 3);
    V.row(0) = e0.transpose();
    V.row(1) = e1.transpose();
    V.row(2) = f0.transpose();
    V.row(3) = f1.transpose();

    Eigen::MatrixX3d X = fix.build_X(V);
    Eigen::MatrixX3d tangents = fix.build_tangents(X);

    const auto [val, grad, hess] = fix.edge3->smooth_edge3_tangent_term_hessian(
        dn, tangents, alpha_t, beta_t);

    Eigen::VectorXd x = Edge3TestFixture::flatten_with_dir(dn, tangents);

    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        x,
        [&](const Eigen::VectorXd& y) {
            Eigen::MatrixX3d T = fd::unflatten(y.tail(y.size() - 3), 3);
            return std::get<0>(fix.edge3->smooth_edge3_tangent_term_gradient(
                y.head(3), T, alpha_t, beta_t));
        },
        fgrad, fd::AccuracyOrder::FOURTH, 1e-5);

    CHECK((fgrad - grad).norm() / 1e-8 <= grad.norm());

    Eigen::MatrixXd fhess;
    fd::finite_jacobian(
        x,
        [&](const Eigen::VectorXd& y) -> Eigen::VectorXd {
            Eigen::MatrixX3d T = fd::unflatten(y.tail(y.size() - 3), 3);
            return std::get<1>(fix.edge3->smooth_edge3_tangent_term_gradient(
                y.head(3), T, alpha_t, beta_t));
        },
        fhess, fd::AccuracyOrder::SECOND, 1e-7);

    CHECK((fhess - hess).norm() / 1e-6 <= hess.norm());
}

// =============================================================================
// Normal term hessian — 2 VARIANT neighbors (coverage for autodiff inner loop)
// =============================================================================
TEST_CASE(
    "Edge normal term 2 VARIANT neighbors", "[distance][edge-edge][gradient]")
{
    // Use parameters that make normal types VARIANT for the test geometry.
    const double alpha_n = 0.85, beta_n = 0.2;
    SmoothContactParameters params(1e-3, 1, 0, alpha_n, beta_n, 2);

    // dn is the raw direction from the edge to the outside point (d).
    // The normal-term functions expect direction = -d.normalized(), so we
    // pass neg_dn = -dn below.
    Eigen::Vector3d e0(0, 0, 0), e1(1, 0, 0), dn(0, 0, 1);
    Eigen::Vector3d neg_dn = -dn; // direction = -d.normalized()

    // Face normals for the fixture winding are:
    //   face 0 (faces=[0,1,2]): n = (0, -fz0, fy0), so dn·n_hat = fy0/|n|
    //   face 1 (faces=[1,0,3]): n = (0, fz1, -fy1), so dn·n_hat = -fy1/|n|
    // To get VARIANT types (dn·n_hat ∈ (-0.85, 0.2)) with normal_sum < 1,
    // use f0 with negative y and f1 with positive y, and large enough z to
    // push the dot products well into the transition region.
    Eigen::Vector3d f0(0.4, -0.3, GENERATE(take(5, random(0.4, 0.8))));
    Eigen::Vector3d f1(0.6, 0.3, GENERATE(take(5, random(0.4, 0.8))));

    Eigen::MatrixX3d fv(2, 3);
    fv.row(0) = f0.transpose();
    fv.row(1) = f1.transpose();

    // Build fixture and mesh
    auto fix = Edge3TestFixture::build(e0, e1, fv, dn, params, true);
    REQUIRE(fix.edge3->n_face_neighbors() == 2);

    Eigen::MatrixXd V(4, 3);
    V.row(0) = e0.transpose();
    V.row(1) = e1.transpose();
    V.row(2) = f0.transpose();
    V.row(3) = f1.transpose();

    Eigen::MatrixX3d X = fix.build_X(V);

    // Flatten [neg_dn, X] into a single vector for FD
    Eigen::VectorXd x = Edge3TestFixture::flatten_with_dir(neg_dn, X);

    // Evaluate analytic hessian via Edge3 (direction = -d.normalized() =
    // neg_dn)
    const auto [val, grad, hess] =
        fix.edge3->smooth_edge3_normal_term_hessian(neg_dn, X, alpha_n, beta_n);

    // Basic sanity: value should be finite and between 0 and 1.
    CHECK(val >= 0.0);
    CHECK(val <= 1.0);

    REQUIRE(grad.isZero() == false);
    REQUIRE(hess.isZero() == false);

    // Finite-difference gradient
    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        x,
        [&](const Eigen::VectorXd& y) {
            Eigen::MatrixX3d Y = fd::unflatten(y.tail(y.size() - 3), 3);
            return std::get<0>(fix.edge3->smooth_edge3_normal_term_gradient(
                y.head(3), Y, alpha_n, beta_n));
        },
        fgrad, fd::AccuracyOrder::FOURTH, 1e-5);

    // Compare gradient (use same style/tolerance as other tests)
    CHECK((fgrad - grad).norm() / 1e-8 <= grad.norm());

    // Finite-difference Hessian (via Jacobian of gradient)
    Eigen::MatrixXd fhess;
    fd::finite_jacobian(
        x,
        [&](const Eigen::VectorXd& y) -> Eigen::VectorXd {
            Eigen::MatrixX3d Y = fd::unflatten(y.tail(y.size() - 3), 3);
            return std::get<1>(fix.edge3->smooth_edge3_normal_term_gradient(
                y.head(3), Y, alpha_n, beta_n));
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