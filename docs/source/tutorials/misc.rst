Miscellaneous
=============

Static Intersection Checks
--------------------------

Static intersection checks are not a core part of the IPC algorithm, but they are useful for debugging and verifying your input is intersection-free. The IPC Toolkit provides a function for checking if a mesh is intersection-free:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            bool is_intersecting = ipc::has_intersections(collision_mesh, vertices);

    .. md-tab-item:: Python

        .. code-block:: python

            is_intersecting = ipctk.has_intersections(collision_mesh, vertices)



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

.. _vertex-derivative-layout:

Vertex Derivative Layout
------------------------

By default, the IPC Toolkit orders derivative vectors (gradients and Hessians) with respect to vertex DOFs in a **row-major** layout:

.. math::

    \mathbf{x} = [x_0, y_0, z_0, x_1, y_1, z_1, \ldots]

Some simulation frameworks instead use a **column-major** layout, where all coordinates of the same dimension are grouped together:

.. math::

    \mathbf{x} = [x_0, x_1, \ldots, y_0, y_1, \ldots, z_0, z_1, \ldots]

The IPC Toolkit provides a compile-time CMake option to select which layout is used for all gradient and Hessian assembly. This is controlled by the ``IPC_TOOLKIT_VERTEX_DERIVATIVE_LAYOUT`` cache variable, which accepts either ``RowMajor`` (default) or ``ColMajor``.

To build with column-major derivative ordering, pass the option when configuring CMake:

.. code-block:: bash

    cmake -DIPC_TOOLKIT_VERTEX_DERIVATIVE_LAYOUT=ColMajor -B build -S .

Or to explicitly use the default row-major layout:

.. code-block:: bash

    cmake -DIPC_TOOLKIT_VERTEX_DERIVATIVE_LAYOUT=RowMajor -B build -S .

.. note::

    This is an advanced option marked as such in CMake. Most users will not need to change it. It is primarily useful when integrating the IPC Toolkit into a simulation framework that stores DOFs in a column-major layout, avoiding the need to permute vectors and matrices at the interface boundary.

This setting affects all functions that assemble local per-element gradients and Hessians into global vectors and matrices, including ``local_gradient_to_global_gradient``, ``local_hessian_to_global_triplets``, and ``local_jacobian_to_global_triplets``. It also affects ``CollisionMesh::vertex_matrix_to_dof_matrix``, which converts vertex-indexed sparse matrices to DOF-indexed sparse matrices.

The chosen layout is stored as the compile-time constant ``ipc::VERTEX_DERIVATIVE_LAYOUT`` (equal to ``Eigen::RowMajor`` or ``Eigen::ColMajor``) in the generated ``config.hpp`` header.

.. warning::

    This feature is experimental. If you encounter any bugs, please report them on `GitHub <https://github.com/ipc-sim/ipc-toolkit/issues>`_.