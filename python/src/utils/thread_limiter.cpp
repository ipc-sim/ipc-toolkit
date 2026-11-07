#include <common.hpp>

#include <memory>
#include <cstdlib>

#include <tbb/info.h>
#include <tbb/global_control.h>

#include <ipc/utils/logger.hpp>

using namespace ipc;

static std::shared_ptr<tbb::global_control> thread_limiter;

int get_num_threads()
{
    return tbb::global_control::active_value(
        tbb::global_control::max_allowed_parallelism);
}

void set_num_threads(int nthreads)
{
    if (nthreads <= 0) {
        nthreads = tbb::info::default_concurrency();
    } else if (nthreads > tbb::info::default_concurrency()) {
        logger().warn(
            "Attempting to use more threads than available ({:d} > {:d})!",
            nthreads, tbb::info::default_concurrency());
        nthreads = tbb::info::default_concurrency();
    }
    thread_limiter = std::make_shared<tbb::global_control>(
        tbb::global_control::max_allowed_parallelism, nthreads);
}

void define_thread_limiter(py::module_& m)
{
    m.def(
        "get_num_threads", &get_num_threads,
        "get maximum number of threads to use");
    m.def(
        "set_num_threads", &set_num_threads,
        "set maximum number of threads to use", "nthreads"_a);

    // Allow the thread limit to be set globally via an environment variable
    // (e.g., in production) rather than requiring every script to call
    // set_num_threads() itself.
    const char* env_val = std::getenv("TBB_NUM_THREADS");
    if (env_val != nullptr) {
        try {
            set_num_threads(std::stoi(env_val));
        } catch (const std::exception& e) {
            logger().error("Invalid value for TBB_NUM_THREADS: {}", env_val);
        }
    }

    std::atexit([]() { thread_limiter.reset(); });
}
