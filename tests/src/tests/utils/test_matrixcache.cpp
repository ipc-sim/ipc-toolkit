#include <iostream>
#include <catch2/catch_test_macros.hpp>
#include <ipc/utils/local_to_global.hpp>
#include <catch2/catch_approx.hpp>

using namespace ipc;

TEST_CASE("SparseMatrixCache", "[utils][matrix_cache]")
{
    const size_t N = 5;

    {
        auto cache = std::make_shared<SparseMatrixCache>(N);
        CHECK(cache);

        cache = std::make_shared<SparseMatrixCache>(N, 2 * N);
        CHECK(cache);

        cache->init(SparseMatrixCache(N, N), true);
        CHECK(cache);

        cache->set_zero();
        CHECK(
            cache->get_matrix().squaredNorm()
            == Catch::Approx(0.).margin(1e-15));

        SparseMatrixCache cache2(*cache);
    }

    {
        SparseMatrixCache cache(N, N);

        for (size_t i = 0; i < N - 1; i++) {
            cache.add_value(0, i, i + 1, 1.);
        }

        CHECK(
            cache.get_matrix().squaredNorm()
            == Catch::Approx(N - 1.).margin(1e-15));
    }
}
