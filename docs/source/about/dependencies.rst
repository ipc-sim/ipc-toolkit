Dependencies
============

.. figure:: /_static/graphviz/dependencies.svg
   :align: center

   Default dependencies of the ``ipc::toolkit`` library. Excludes CUDA and Python bindings.

The IPC Toolkit depends on a handful of third-party libraries, which are used to provide various functionality.

**All required dependencies are downloaded through CMake** depending on the build options, and are built automatically when you build the IPC Toolkit. You do not need to install them separately.

These dependencies are all `permissively licensed <license.html>`_, and we list them here to give due credit!

.. list-table::
    :header-rows: 1
    :widths: 15 35 10 40

    * - Name
      - Purpose
      - License
      - Link
    * - Eigen
      - Linear algebra
      - MPL-2.0
      - `eigen.tuxfamily.org <http://eigen.tuxfamily.org/>`_
    * - libigl
      - Geometry functions and predicates
      - MPL-2.0
      - `libigl.github.io <https://libigl.github.io/>`_
    * - oneTBB
      - Multithreading
      - Apache-2.0
      - `github.com/oneapi-src/oneTBB <https://github.com/oneapi-src/oneTBB>`_
    * - Tight-Inclusion
      - Provably conservative CCD
      - MIT
      - `github.com/Continuous-Collision-Detection/Tight-Inclusion <https://github.com/Continuous-Collision-Detection/Tight-Inclusion>`_
    * - Scalable-CCD
      - Scalable (GPU) CCD
      - Apache-2.0
      - `github.com/Continuous-Collision-Detection/Scalable-CCD <https://github.com/Continuous-Collision-Detection/Scalable-CCD>`_
    * - spdlog
      - Logger
      - MIT
      - `github.com/gabime/spdlog <https://github.com/gabime/spdlog>`_
    * - TinyAD
      - Automatic differentiation for testing and in non-performance critical smooth contact functions
      - MIT
      - `github.com/microsoft/TinyAD <https://github.com/patr-schm/TinyAD>`_

Optional Dependencies
---------------------

Additionally, IPC Toolkit may optionally use the following libraries:

.. list-table::
    :header-rows: 1

    * - Name
      - Purpose
      - License
      - Link
      - Enabled
      - CMake Option
    * - xsimd
      - Cross-platform SIMD library for vectorization
      - BSD-3-Clause
      - `xsimd.readthedocs.io <https://xsimd.readthedocs.io/en/latest>`_
      - |:white_check_mark:|
      - ``IPC_TOOLKIT_WITH_SIMD``
    * - robin-map
      - Faster hashing
      - MIT
      - `github.com/Tessil/robin-map <https://github.com/Tessil/robin-map>`_
      - |:white_check_mark:|
      - ``IPC_TOOLKIT_WITH_ROBIN_MAP``
    * - Abseil
      - Hashing utilities
      - Apache-2.0
      - `abseil.io <https://abseil.io/>`_
      - |:white_check_mark:|
      - ``IPC_TOOLKIT_WITH_ABSEIL``
    * - filib
      - Interval arithmetic for nonlinear trajectories/CCD
      - LGPL-2.1
      - `github.com/zfergus/filib <https://github.com/zfergus/filib>`_
      - |:white_check_mark:|
      - ``IPC_TOOLKIT_WITH_FILIB``
    * - nlohmann/json
      - JSON parsing for profiler and tests
      - MIT
      - `github.com/nlohmann/json <https://github.com/nlohmann/json>`_
      - |:white_large_square:|
      - ``IPC_TOOLKIT_WITH_PROFILER``
    * - rational-cpp
      - Rational arithmetic used for exact intersection checks (requires `GMP <https://gmplib.org>`_ to be installed at a system level)
      - MIT
      - `github.io/zfergus/rational-cpp <https://github.io/zfergus/rational-cpp>`_
      - |:white_large_square:|
      - ``IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION``
    * - Etienne Vouga's Collision Detection Library
      - Inexact CCD (included for comparison with the original IPC library)
      - ???
      - `github.com/evouga/collisiondetection <https://github.com/evouga/collisiondetection>`_
      - |:white_large_square:|
      - ``IPC_TOOLKIT_WITH_INEXACT_CCD``

Some of these libraries are enabled by default, and some are not. You can enable or disable them by passing the appropriate CMake option when you configure the IPC Toolkit build.

.. warning::
    ``filib`` is licensed under `LGPL-2.1 <https://github.com/zfergus/filib/blob/main/LICENSE>`_ and as such it is required to be dynamically linked. Doing so automatically is a challenge, so by default we use static linkage. Enabling dynamic linkage requires copying the ``.so``/``.dylib``/``.dll`` file to the binary directory or system path. To enable this, set the CMake option ``FILIB_BUILD_SHARED_LIBS`` to ``ON`` and add this CMake code to copy the shared library object to the binary directory:

    .. code-block:: cmake

        # Copy shared lib to the output directory
        add_custom_command(
            TARGET ${MY_EXE_TARGET} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:${MY_EXE_TARGET}> $<TARGET_FILE_DIR:${MY_EXE_TARGET}>
            COMMAND_EXPAND_LISTS
        )

    where ``${MY_EXE_TARGET}`` is the name of your executable target. If you know a better way to handle this, please `let us know <https://github.com/ipc-sim/ipc-toolkit/discussions>`_!

    If you would rather avoid LGPL code entirely, you can disable filib by setting ``IPC_TOOLKIT_WITH_FILIB`` to ``OFF``. With this option disabled, CMake will not download or use any of filib's code.

Python Bindings Dependencies
----------------------------

The following dependencies are optionally used when building Python bindings:

.. list-table::
    :header-rows: 1

    * - Name
      - Purpose
      - License
      - Link
      - Enabled
      - CMake Option
    * - Pybind11
      - Python bindings
      - BSD-3-Clause
      - `github.com/pybind/pybind11 <https://github.com/pybind/pybind11>`_
      - |:white_check_mark:|
      - ``IPC_TOOLKIT_BUILD_PYTHON``
    * - pybind11_json
      - JSON support for pybind11 (used in profiler)
      - MIT
      - `github.com/pybind/pybind11_json <https://github.com/pybind/pybind11_json>`_
      - |:white_large_square:|
      - ``IPC_TOOLKIT_BUILD_PYTHON`` and ``IPC_TOOLKIT_WITH_PROFILER``

Unit Test Dependencies
----------------------

The following dependencies are optionally used when unit tests are enabled:

.. list-table::
    :header-rows: 1

    * - Name
      - Purpose
      - License
      - Link
      - Enabled
      - CMake Option
    * - Catch2
      - Testing framework
      - BSL-1.0
      - `github.com/catchorg/Catch2 <https://github.com/catchorg/Catch2.git>`_
      - |:white_check_mark:|
      - ``IPC_TOOLKIT_BUILD_TESTS``
    * - finite-diff
      - Finite-difference comparisons
      - MIT
      - `github.com/zfergus/finite-diff <https://github.com/zfergus/finite-diff>`_
      - |:white_check_mark:|
      - ``IPC_TOOLKIT_BUILD_TESTS``