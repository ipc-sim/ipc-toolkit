#include <pybind11/pybind11.h>

#include <tbb/global_control.h>
#include <tbb/task_scheduler_init.h>

namespace py = pybind11;

// static tbb::global_control thread_limiter(
//     tbb::global_control::max_allowed_parallelism,
//     tbb::task_scheduler_init::default_num_threads());
//
// void set_num_threads(int nthreads)
// {
//     if (nthreads <= 0) {
//         nthreads = tbb::task_scheduler_init::default_num_threads();
//     } else if (nthreads > tbb::task_scheduler_init::default_num_threads()) {
//         logger().warn(
//             "Attempting to use more threads than available ({:d} > "
//             "{:d})!",
//             nthreads, tbb::task_scheduler_init::default_num_threads());
//         nthreads = tbb::task_scheduler_init::default_num_threads();
//     }
//     thread_limiter = tbb::global_control(
//         tbb::global_control::max_allowed_parallelism, nthreads);
// }

void define_thread_limiter_functions(py::module_& m)
{
    // m.def(
    //     "set_num_threads", &set_num_threads, "maximum number of threads to
    //     use", py::arg("nthreads"));
}
