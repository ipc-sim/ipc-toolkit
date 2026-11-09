#include <catch2/catch_test_macros.hpp>

#include <ipc/utils/eigen_abi.hpp>

using namespace ipc;

TEST_CASE(
    "Library and its consumers agree on Eigen's ABI", "[utils][eigen_abi]")
{
    // Code that links the library compiles with the same Eigen settings, so a
    // mismatch here means Eigen objects are unsafe to pass between the two.
    CHECK(eigen_abi_matches_library());
}
