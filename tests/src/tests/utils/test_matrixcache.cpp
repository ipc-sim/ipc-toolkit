#include <iostream>
#include <catch2/catch_test_macros.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <iostream>

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

        for (size_t i = 0; i < N - 1; i++) {
            cache->add_value(0, i, i + 1, 1.);
        }

        CHECK(fabs(cache->get_matrix().squaredNorm() - (N - 1.)) < 1.e-12);
    }

    {
        SparseMatrixCache cache1(N, N);
        SparseMatrixCache cache2(N, N);
        for (size_t i = 0; i < N - 1; i++) {
            cache1.add_value(0, i, i + 1, 1.);
            cache2.add_value(0, i + 1, i, 1.);
        }

        cache2.prune();
        cache1 += cache2;

        CHECK(cache1.get_matrix().squaredNorm() - 2 * (N - 1.) < 1.e-12);
    }
}

TEST_CASE("DenseMatrixCache", "[utils][matrix_cache]")
{
    const size_t N = 5;

    {
        auto cache = std::make_shared<DenseMatrixCache>(N);
        CHECK(cache);

        cache = std::make_shared<DenseMatrixCache>(N, 2 * N);
        CHECK(cache);

        cache->init(DenseMatrixCache(N, N));
        CHECK(cache);

        cache->set_zero();

        for (size_t i = 0; i < N - 1; i++) {
            cache->add_value(0, i, i + 1, 1.);
        }

        CHECK(fabs(cache->get_matrix().squaredNorm() - (N - 1.)) < 1.e-12);
    }

    {
        DenseMatrixCache cache1(N, N);
        DenseMatrixCache cache2(N, N);
        for (size_t i = 0; i < N - 1; i++) {
            cache1.add_value(0, i, i + 1, 1.);
            cache2.add_value(0, i + 1, i, 1.);
        }

        cache2.prune();
        cache1 += cache2;

        CHECK(cache1.get_matrix().squaredNorm() - 2 * (N - 1.) < 1.e-12);
    }
}
