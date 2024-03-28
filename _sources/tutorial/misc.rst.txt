Miscellaneous
=============

Static Intersection Checks
--------------------------

Static intersection checks are not a core part of the IPC algorithm, but they are useful for debugging and verifying your input is intersection-free. The IPC Toolkit provides a function for checking if a mesh is intersection-free:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            bool is_intersecting = ipc::has_intersections(collision_mesh, displaced);

    .. md-tab-item:: Python

        .. code-block:: python

            is_intersecting = ipctk.has_intersections(collision_mesh, displaced)



Logger
------

The IPC Toolkit uses the `spdlog <https://github.com/gabime/spdlog>`_ library for logging. By default, the IPC Toolkit will log to ``stdout``. To change this behavior, you can set the logger using the ``ipc::set_logger`` function. For example, to additionally log to a file, you can do the following:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            std::string log_file = "log.txt";

            std::vector<spdlog::sink_ptr> sinks;
            sinks.emplace_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
            sinks.emplace_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file, /*truncate=*/true));

            ipc::set_logger(std::make_shared<spdlog::logger>("ipctk", sinks.begin(), sinks.end()));

    .. md-tab-item:: Python

        .. code-block:: python

            # Unfortunately, this is not yet supported in Python.

You can also set the log level:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            ipc::logger().set_level(spdlog::level::debug);

    .. md-tab-item:: Python

        .. code-block:: python

            ipctk.set_logger_level(ipctk.LoggerLevel.debug)

Multi-threading
---------------

The IPC Toolkit utilizes `oneTBB <https://oneapi-src.github.io/oneTBB>`_ to enable multi-threading functionality. This is enabled by default, and significantly improves performance.
However, with multi-threading enabled, it is expected that the results can be non-deterministic because of rounding differences when adding numbers in different orders. To make the results deterministic you can limit TBB's maximum number of threads to one. The following code shows how to do this:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            #include <tbb/global_control.h>

            // ...

            tbb::global_control thread_limiter(tbb::global_control::max_allowed_parallelism, 1);

        As long as this `thread_limiter` object stays alive, the number of threads will be limited to 1. This example shows the thread limiter stack allocated, so when the object goes out-of-scope, the limit will be released. If you want this limit for the entire length of your program it is best to define it in ``main`` or to define it statically using a `std::shared_ptr<tbb::global_control>`. For example,

        .. code-block:: c++

            #include <tbb/info.h>
            #include <tbb/global_control.h>

            static std::shared_ptr<tbb::global_control> thread_limiter;

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

    .. md-tab-item:: Python

        .. code-block:: python

            ipctk.set_num_threads(1)

        This limit will persist for the duration of the program or until you call

        .. code-block:: python

            ipctk.set_num_threads(-1)

        to reset the number of threads to the default.

You can also get the current maximum number of threads as follows:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            int nthreads = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

        Additionally, you can get the number of threads TBB will use by default:

        .. code-block:: c++

            #include <tbb/info.h>

            // ...

            int max_nthreads = tbb::info::default_concurrency();

    .. md-tab-item:: Python

        .. code-block:: python

            nthreads = ipctk.get_num_threads()