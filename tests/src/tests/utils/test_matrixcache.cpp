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

    {
        SparseMatrixCache cache1(N, N);
        std::unique_ptr<MatrixCache> cache2 =
            std::make_unique<SparseMatrixCache>(N, N);
        for (size_t i = 0; i < N - 1; i++) {
            cache1.add_value(0, i, i + 1, 1.);
            cache2->add_value(0, i + 1, i, 1.);
        }

        cache2->prune();
        cache1 += *cache2;

        CHECK(
            cache1.get_matrix().squaredNorm()
            == Catch::Approx(2 * (N - 1.)).margin(1e-15));
    }

    {
        SparseMatrixCache cache1(N, N);
        SparseMatrixCache cache2(N, N);
        SparseMatrixCache cache3(N, N);
        for (size_t i = 0; i < N - 1; i++) {
            cache1.add_value(i / 2, i, i + 1, 1.);
            cache2.add_value(i / 2, i + 1, i, 2.);
            cache3.add_value(i / 2, i, i, 1.);
        }

        CHECK(
            cache2.get_matrix().squaredNorm()
            == Catch::Approx(4 * (N - 1.)).margin(1e-15));

        for (size_t i = 0; i < N - 1; i++) {
            cache2.add_value(i / 2, i + 1, i, 1.);
        }
        CHECK(
            cache2.get_matrix().squaredNorm()
            == Catch::Approx(N - 1.).margin(1e-15));

        cache2.prune();
        cache1 += cache2;

        CHECK(
            cache1.get_matrix().squaredNorm()
            == Catch::Approx(2 * (N - 1.)).margin(1e-15));

        cache3.prune();
        cache3 += cache1;

        CHECK(
            cache3.get_matrix().squaredNorm()
            == Catch::Approx(3 * (N - 1.)).margin(1e-15));
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

        CHECK(
            cache->get_matrix().squaredNorm()
            == Catch::Approx(N - 1).margin(1e-15));
    }

    {
        DenseMatrixCache cache1(N, N);
        std::unique_ptr<MatrixCache> cache2 =
            std::make_unique<DenseMatrixCache>(N, N);
        for (size_t i = 0; i < N - 1; i++) {
            cache1.add_value(0, i, i + 1, 1.);
            cache2->add_value(0, i + 1, i, 1.);
        }

        cache2->prune();
        cache1 += *cache2;

        CHECK(
            cache1.get_matrix().squaredNorm()
            == Catch::Approx(2 * (N - 1.)).margin(1e-15));
    }
}
