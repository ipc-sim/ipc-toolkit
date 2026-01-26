#include <ipc/utils/profiler.hpp>

#ifdef IPC_TOOLKIT_WITH_PROFILER

#include <catch2/catch_all.hpp>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

using namespace ipc;

TEST_CASE("Profiler", "[profiler]")
{
    constexpr int sleep_time_ms = 100;
    constexpr int num_threads = 10;

    {
        IPC_TOOLKIT_PROFILE_BLOCK("Block 1");
        std::this_thread::sleep_for(std::chrono::milliseconds(sleep_time_ms));
    }

    nlohmann::json data = profiler().data();

    CHECK(data.size() == 1);

    REQUIRE(data.contains("Block 1"));
    nlohmann::json block1 = data.at("Block 1");

    CHECK(block1.size() == 2); // count, time_ms
    CHECK(
        block1["time_ms"].get<double>()
        == Catch::Approx(sleep_time_ms).margin(10));
    CHECK(block1["count"].get<int>() == 1);

    // ---------------------------------------------------------------------

    tbb::parallel_for(0, num_threads, [&](int) {
        // for (int i = 0; i < num_threads; ++i) {
        IPC_TOOLKIT_PROFILE_BLOCK("Block 2");
        std::this_thread::sleep_for(
            std::chrono::milliseconds(sleep_time_ms / num_threads));
        // }
    });

    data = profiler().data();
    profiler().print();

    CHECK(data.size() == 2);

    REQUIRE(data.contains("Block 2"));
    nlohmann::json block2 = data.at("Block 2");
    CHECK(block2.size() == 2); // count, time_ms
    CHECK(block2["count"].get<int>() == num_threads);

    // ---------------------------------------------------------------------

    {
        IPC_TOOLKIT_PROFILE_BLOCK("Block 3");

        tbb::parallel_for(0, num_threads, [&](int) {
            // for (int i = 0; i < num_threads; ++i) {
            IPC_TOOLKIT_PROFILE_BLOCK("Block 4");
            {
                IPC_TOOLKIT_PROFILE_BLOCK("Block 5");
                std::this_thread::sleep_for(
                    std::chrono::milliseconds(sleep_time_ms / num_threads));
            }
            // }
        });
    }

    data = profiler().data();
    profiler().print();

    CHECK(data.size() == 3);

    REQUIRE(data.contains("Block 3"));
    CAPTURE(!data.contains("Block 4"));
    nlohmann::json block3 = data.at("Block 3");
    CHECK(block3.size() == 3); // count, time_ms, Block 4
    CHECK(block3["count"].get<int>() == 1);

    REQUIRE(block3.contains("Block 4"));
    nlohmann::json block4 = block3.at("Block 4");
    CHECK(block4.size() == 3); // count, time_ms
    CHECK(block4["count"].get<int>() == num_threads);
}

#endif