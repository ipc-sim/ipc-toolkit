#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <ipc/distance/point_point.hpp>
#include <ipc/distance/line_line.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/barrier/barrier.hpp>
#include <ipc/distance/edge_edge.hpp>

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

// 1. The fast worker (Compiler knows N)
template <int N, typename DerivedA, typename DerivedB>
inline double fast_dist(
    const Eigen::MatrixBase<DerivedA>& p0,
    const Eigen::MatrixBase<DerivedB>& p1)
{
    return (p1.template head<N>() - p0.template head<N>()).squaredNorm();
}

// 2. The dispatcher
double point_point_distance_dispatch(
    Eigen::ConstRef<VectorMax3d> p0, Eigen::ConstRef<VectorMax3d> p1)
{
    // By the time we reach the worker, N is a constant
    if (p0.size() == 3) {
        return fast_dist<3>(p0, p1);
    }
    if (p0.size() == 2) {
        return fast_dist<2>(p0, p1);
    }
    return (p1 - p0).squaredNorm();
}

// 2. The dispatcher
template <int N, int M = N>
double point_point_distance_templated(
    Eigen::Matrix<double, N, 1, Eigen::ColMajor, M, 1> p0,
    Eigen::Matrix<double, N, 1, Eigen::ColMajor, M, 1> p1)
{
    return (p1 - p0).squaredNorm();
}

inline double point_point_distance_templated(
    Eigen::ConstRef<VectorMax3d> p0, Eigen::ConstRef<VectorMax3d> p1)
{
    if (p0.size() == 3) {
        return point_point_distance_templated<3, 3>(p0, p1);
    }
    // else if (p0.size() == 2)
    return point_point_distance_templated<2, 2>(p0, p1);
    // return point_point_distance_templated<-1, 3>(p0, p1);
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

    BENCHMARK("Dispatch Nd")
    {
        return point_point_distance_dispatch(V.row(p0i), V.row(p1i));
    };

    BENCHMARK("Templated")
    {
        return point_point_distance_templated<3, 3>(V.row(p0i), V.row(p1i));
    };

    BENCHMARK("Templated(Nd)")
    {
        return point_point_distance_templated(V.row(p0i), V.row(p1i));
    };
}

// =============================================================================

double line_line_distance_ref(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    const Eigen::Vector3d normal = (ea1 - ea0).cross(eb1 - eb0);
    const double line_to_line = (eb0 - ea0).dot(normal);
    return line_to_line * line_to_line / normal.squaredNorm();
}

template <
    typename DerivedA,
    typename DerivedB,
    typename DerivedC,
    typename DerivedD>
double line_line_distance_templated(
    const Eigen::MatrixBase<DerivedA>& ea0,
    const Eigen::MatrixBase<DerivedB>& ea1,
    const Eigen::MatrixBase<DerivedC>& eb0,
    const Eigen::MatrixBase<DerivedD>& eb1)
{
    const Eigen::Vector3d normal = (ea1 - ea0).template block<3, 1>(0, 0).cross(
        (eb1 - eb0).template block<3, 1>(0, 0));
    const double line_to_line = (eb0 - ea0).dot(normal);
    return line_to_line * line_to_line / normal.squaredNorm();
}

double line_line_distance_row_ref(
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
    const int N = 100;
    const Eigen::MatrixXd V = Eigen::MatrixXd::Random(N, 3);
    const Eigen::Matrix<double, -1, -1, Eigen::RowMajor> VR = V;

    // Use a heap-allocated vector of random indices to defeat LICM. The prime
    // size means i % M has no power-of-2 period, and the heap pointer prevents
    // the compiler from tracing provenance back to V. Catch2 has no
    // DoNotOptimize, so this is the most portable approach.
    constexpr int M = 997; // prime
    Eigen::ArrayXi idx =
        Eigen::ArrayXi::Random(M).abs().unaryExpr([&](int x) { return x % N; });

    BENCHMARK("Templated", i)
    {
        return line_line_distance_templated(
            V.row(idx[i % M]), V.row(idx[(i + 1) % M]), V.row(idx[(i + 2) % M]),
            V.row(idx[(i + 3) % M]));
    };

    BENCHMARK("Lib", i)
    {
        return line_line_distance(
            V.row(idx[i % M]), V.row(idx[(i + 1) % M]), V.row(idx[(i + 2) % M]),
            V.row(idx[(i + 3) % M]));
    };

    BENCHMARK("Ref", i)
    {
        return line_line_distance_ref(
            V.row(idx[i % M]), V.row(idx[(i + 1) % M]), V.row(idx[(i + 2) % M]),
            V.row(idx[(i + 3) % M]));
    };

    BENCHMARK("Ref(Copy)", i)
    {
        const Eigen::Vector3d ea0 = V.row(idx[(i + 0) % M]);
        const Eigen::Vector3d ea1 = V.row(idx[(i + 1) % M]);
        const Eigen::Vector3d eb0 = V.row(idx[(i + 2) % M]);
        const Eigen::Vector3d eb1 = V.row(idx[(i + 3) % M]);
        return line_line_distance_ref(ea0, ea1, eb0, eb1);
    };

    BENCHMARK("Ref (VR)", i)
    {
        return line_line_distance(
            VR.row(idx[i % M]), VR.row(idx[(i + 1) % M]),
            VR.row(idx[(i + 2) % M]), VR.row(idx[(i + 3) % M]));
    };

    BENCHMARK("RowRef", i)
    {
        return line_line_distance_row_ref(
            V.row(idx[i % M]), V.row(idx[(i + 1) % M]), V.row(idx[(i + 2) % M]),
            V.row(idx[(i + 3) % M]));
    };

    BENCHMARK("Copy", i)
    {
        return line_line_distance_copy(
            V.row(idx[i % M]), V.row(idx[(i + 1) % M]), V.row(idx[(i + 2) % M]),
            V.row(idx[(i + 3) % M]));
    };
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

// =============================================================================

TEST_CASE("Float vs Double: Distance", "[!benchmark][eigen][float]")
{
    constexpr int N = 100, M = 1000;
    const Eigen::MatrixXf Vf = Eigen::MatrixXf::Random(N, 3);
    const Eigen::MatrixXd Vd = Vf.cast<double>();

    BENCHMARK("float")
    {
        float sum = 0.0;
        for (int i = 0; i < M; ++i) {
            sum += barrier(
                edge_edge_distance<float>(
                    Vf.row(i % N), Vf.row((50 + i) % N), Vf.row((75 + i) % N),
                    Vf.row((99 + i) % N)),
                1.0f);
        }
        return sum;
    };
    BENCHMARK("double")
    {
        double sum = 0.0;
        for (int i = 0; i < M; ++i) {
            sum += barrier(
                edge_edge_distance<double>(
                    Vd.row(i % N), Vd.row((50 + i) % N), Vd.row((75 + i) % N),
                    Vd.row((99 + i) % N)),
                1.0);
        }
        return sum;
    };
}

TEST_CASE("Float vs Double: Hessian", "[!benchmark][eigen][float]")
{
    constexpr int N = 100, M = 1000;
    const Eigen::MatrixXf Vf = Eigen::MatrixXf::Random(N, 3);
    const Eigen::MatrixXd Vd = Vf.cast<double>();

    BENCHMARK("float")
    {
        Eigen::Matrix<float, 12, 12> hess =
            Eigen::Matrix<float, 12, 12>::Zero();
        for (int i = 0; i < M; ++i) {
            auto d = edge_edge_distance<float>(
                Vf.row(i % N), Vf.row((50 + i) % N), Vf.row((75 + i) % N),
                Vf.row((99 + i) % N));
            auto grad_d = edge_edge_distance_gradient<float>(
                Vf.row(i % N), Vf.row((50 + i) % N), Vf.row((75 + i) % N),
                Vf.row((99 + i) % N));
            auto hess_d = edge_edge_distance_hessian<float>(
                Vf.row(i % N), Vf.row((50 + i) % N), Vf.row((75 + i) % N),
                Vf.row((99 + i) % N));

            auto db = barrier_first_derivative<float>(d, 1.0f);
            auto d2b = barrier_second_derivative<float>(d, 1.0f);

            hess += db * hess_d + (d2b * grad_d) * grad_d.transpose();
        }
        return hess;
    };
    BENCHMARK("double")
    {
        Eigen::Matrix<double, 12, 12> hess =
            Eigen::Matrix<double, 12, 12>::Zero();
        for (int i = 0; i < M; ++i) {
            auto d = edge_edge_distance<double>(
                Vd.row(i % N), Vd.row((50 + i) % N), Vd.row((75 + i) % N),
                Vd.row((99 + i) % N));
            auto grad_d = edge_edge_distance_gradient<double>(
                Vd.row(i % N), Vd.row((50 + i) % N), Vd.row((75 + i) % N),
                Vd.row((99 + i) % N));
            auto hess_d = edge_edge_distance_hessian<double>(
                Vd.row(i % N), Vd.row((50 + i) % N), Vd.row((75 + i) % N),
                Vd.row((99 + i) % N));

            auto db = barrier_first_derivative<double>(d, 1.0);
            auto d2b = barrier_second_derivative<double>(d, 1.0);

            hess += db * hess_d + (d2b * grad_d) * grad_d.transpose();
        }
        return hess;
    };
}
