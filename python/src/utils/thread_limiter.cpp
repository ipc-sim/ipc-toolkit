#include <pybind11/pybind11.h>

#include <memory>

#include <tbb/info.h>
#include <tbb/global_control.h>

#include <ipc/utils/logger.hpp>

namespace py = pybind11;

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
        IPC_LOG(warn(
            "Attempting to use more threads than available ({:d} > {:d})!",
            nthreads, tbb::info::default_concurrency()));
        nthreads = tbb::info::default_concurrency();
    }
    thread_limiter = std::make_shared<tbb::global_control>(
        tbb::global_control::max_allowed_parallelism, nthreads);
}

void define_thread_limiter_functions(py::module_& m)
{
    m.def(
        "get_num_threads", &get_num_threads,
        "get maximum number of threads to use");
    m.def(
        "set_num_threads", &set_num_threads,
        "set maximum number of threads to use", py::arg("nthreads"));
}
