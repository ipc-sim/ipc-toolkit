Miscellaneous
=============

Static Intersection Checks
--------------------------

Static intersection checks are not a core part of the IPC algorithm, but they are useful for debugging and verifying your input is intersection free. The IPC Toolkit provides a function for checking if a mesh is intersection free:

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

            # Unfourtnately, this is not yet supported in Python.

You can also set the log level:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            ipc::logger().set_level(spdlog::level::debug);

    .. md-tab-item:: Python

        .. code-block:: python

            ipctk.set_logger_level(ipctk.LoggerLevel.debug)
