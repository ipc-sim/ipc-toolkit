#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <ipc/distance/point_point.hpp>
#include <ipc/distance/line_line.hpp>
#include <ipc/distance/point_triangle.hpp>

#include <Eigen/Geometry>

using namespace ipc;

TEST_CASE("Template dynamic vs static", "[!benchmark][eigen]")
{
    const Eigen::MatrixXd V = Eigen::MatrixXd::Random(100, 3);

    int vi = 0, t0i = 50, t1i = 75, t2i = 99;

    BENCHMARK("Static")
    {
        Eigen::Vector3d p = V.row(vi);
        Eigen::Vector3d t0 = V.row(t0i);
        Eigen::Vector3d t1 = V.row(t1i);
        Eigen::Vector3d t2 = V.row(t2i);

        return point_triangle_distance_hessian(p, t0, t1, t2);
    };

    BENCHMARK("Dynamic")
    {
        return point_triangle_distance_hessian(
            V.row(vi), V.row(t0i), V.row(t1i), V.row(t2i));
    };
}

// =============================================================================

double point_point_distance_ref(
    const Eigen::Ref<const Eigen::Vector3d>& p0,
    const Eigen::Ref<const Eigen::Vector3d>& p1)
{
    return (p1 - p0).squaredNorm();
}

double
point_point_distance_copy(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1)
{
    return (p1 - p0).squaredNorm();
}

double
point_point_distance_copy_Nd(const VectorMax3d& p0, const VectorMax3d& p1)
{
    return (p1 - p0).squaredNorm();
}

TEST_CASE("Templated vs Eigen::Ref (PP)", "[!benchmark][eigen][pp][ref]")
{
    Eigen::MatrixXd V = Eigen::MatrixXd::Random(100, 3);

    const int p0i = 0, p1i = V.rows() - 1;

    BENCHMARK("Ref")
    {
        return point_point_distance_ref(V.row(p0i), V.row(p1i));
    };

    BENCHMARK("Ref Copied 3D")
    {
        const Eigen::Vector3d p0 = V.row(p0i);
        const Eigen::Vector3d p1 = V.row(p1i);
        return point_point_distance_ref(p0, p1);
    };

    BENCHMARK("Ref Copied ND")
    {
        const VectorMax3d p0_ND = V.row(p0i);
        const VectorMax3d p1_ND = V.row(p1i);
        return point_point_distance_ref(p0_ND, p1_ND);
    };

    BENCHMARK("Ref Nd")
    {
        return point_point_distance(V.row(p0i), V.row(p1i));
    };

    BENCHMARK("Copy")
    {
        return point_point_distance_copy(V.row(p0i), V.row(p1i));
    };

    BENCHMARK("Copy Nd")
    {
        return point_point_distance_copy_Nd(V.row(p0i), V.row(p1i));
    };
}

// =============================================================================

double line_line_distance_ref(
    const Eigen::Ref<const Eigen::Vector3d, 0, Eigen::InnerStride<>>& ea0,
    const Eigen::Ref<const Eigen::Vector3d, 0, Eigen::InnerStride<>>& ea1,
    const Eigen::Ref<const Eigen::Vector3d, 0, Eigen::InnerStride<>>& eb0,
    const Eigen::Ref<const Eigen::Vector3d, 0, Eigen::InnerStride<>>& eb1)
{
    const Eigen::Vector3d normal = (ea1 - ea0).cross(eb1 - eb0);
    const double line_to_line = (eb0 - ea0).dot(normal);
    return line_to_line * line_to_line / normal.squaredNorm();
}

double line_line_distance_copy(
    const Eigen::Vector3d& ea0,
    const Eigen::Vector3d& ea1,
    const Eigen::Vector3d& eb0,
    const Eigen::Vector3d& eb1)
{
    const Eigen::Vector3d normal = (ea1 - ea0).cross(eb1 - eb0);
    const double line_to_line = (eb0 - ea0).dot(normal);
    return line_to_line * line_to_line / normal.squaredNorm();
}

TEST_CASE("Templated vs Eigen::Ref (LL)", "[!benchmark][eigen][ll][ref]")
{
    const Eigen::MatrixXd V = Eigen::MatrixXd::Random(100, 3);

    int ea0i = 0, ea1i = 50, eb0i = 75, eb1i = 99;

    const auto ea0 = V.row(ea0i);
    const auto ea1 = V.row(ea1i);
    const auto eb0 = V.row(eb0i);
    const auto eb1 = V.row(eb1i);

    BENCHMARK("Templated") { return line_line_distance(ea0, ea1, eb0, eb1); };

    BENCHMARK("Ref") { return line_line_distance_ref(ea0, ea1, eb0, eb1); };

    BENCHMARK("Copy") { return line_line_distance_copy(ea0, ea1, eb0, eb1); };
}

// =============================================================================

template <typename DerivedHess>
void line_line_distance_hessian(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1,
    Eigen::PlainObjectBase<DerivedHess>& hess)
{
    hess.resize(
        ea0.size() + ea1.size() + eb0.size() + eb1.size(),
        ea0.size() + ea1.size() + eb0.size() + eb1.size());
    autogen::line_line_distance_hessian(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], hess.data());
}

TEST_CASE("Return type", "[!benchmark][eigen][return]")
{
    const Eigen::MatrixXd V = Eigen::MatrixXd::Random(100, 3);

    int ea0i = 0, ea1i = 50, eb0i = 75, eb1i = 99;

    const Eigen::Vector3d ea0 = V.row(ea0i);
    const Eigen::Vector3d ea1 = V.row(ea1i);
    const Eigen::Vector3d eb0 = V.row(eb0i);
    const Eigen::Vector3d eb1 = V.row(eb1i);

    BENCHMARK("Explicit")
    {
        const Matrix12d hess = line_line_distance_hessian(ea0, ea1, eb0, eb1);
        return hess;
    };

    BENCHMARK("Implicit")
    {
        Matrix12d hess;
        line_line_distance_hessian(ea0, ea1, eb0, eb1, hess);
        return hess;
    };
}
