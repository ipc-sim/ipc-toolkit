#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>

#include <ipc/utils/math.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
#include <iostream>

Eigen::Matrix3i fun1()
{
    Eigen::Matrix3i A;
    A << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    return A;
}

std::tuple<int, Eigen::Matrix3i> func2()
{
    Eigen::Matrix3i A;
    A << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    return std::make_tuple(0., A);
}

TEST_CASE("assign sub matrix", "[eigen_unit]")
{
    Eigen::Matrix4i B;
    B.setZero();
    B.bottomRightCorner<3, 3>() = fun1();

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            CHECK(B(i + 1, j + 1) == 3 * i + j + 1);
}

TEST_CASE("assign sub matrix in tuple", "[eigen_unit]")
{
    double b;
    Eigen::Matrix4i B;
    B.setZero();
    Eigen::Ref<Eigen::Matrix3i> ref = B.bottomRightCorner<3, 3>();
    std::tie(b, ref) = func2();

    CHECK(b == 0);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            CHECK(B(i + 1, j + 1) == 3 * i + j + 1);
}

TEST_CASE("eigen map and ref", "[eigen_unit]")
{
    Eigen::Matrix4i B;
    B.setZero();
    Eigen::Ref<Eigen::Vector<int, 12>> ref =
        Eigen::Map<Eigen::Vector<int, 12>>(B.data() + 4, 12);
    ref(5) = 1;

    std::cout << "B:\n" << B << "\n";
    std::cout << "ref:\n" << ref.transpose() << "\n";
}

TEST_CASE("eigen ref and seqN", "[eigen_unit]")
{
    Eigen::Matrix4i B;
    B.setZero();
    B({ 1, 2, 3 }, { 1, 2, 3 }) = fun1();

    std::cout << "B:\n" << B << "\n";
}
