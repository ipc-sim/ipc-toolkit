C++
===

.. role:: cpp(code)
   :language: c++
.. role:: cmake(code)
   :language: cmake

.. image:: https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml/badge.svg
   :target: https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml
   :alt: Build
.. image:: https://github.com/ipc-sim/ipc-toolkit/actions/workflows/docs.yml/badge.svg
   :target: https://ipctk.xyz/
   :alt: Docs
.. image:: https://codecov.io/github/ipc-sim/ipc-toolkit/graph/badge.svg?token=9BR6GPKRY8
   :target: https://codecov.io/github/ipc-sim/ipc-toolkit
.. image:: https://img.shields.io/github/license/ipc-sim/ipc-toolkit.svg?color=blue
   :target: https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE
   :alt: License

Build
-----

The easiest way to add the toolkit to an existing CMake project is to download it through CMake.
CMake provides functionality for doing this called `FetchContent <https://cmake.org/cmake/help/latest/module/FetchContent.html>`__ (requires CMake ≥ 3.14).
We use this same process to download all external dependencies.

For example,

.. code:: cmake

   include(FetchContent)
   FetchContent_Declare(
       ipc_toolkit
       GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
       GIT_TAG ${IPC_TOOLKIT_GIT_TAG}
   )
   FetchContent_MakeAvailable(ipc_toolkit)

where :cmake:`IPC_TOOLKIT_GIT_TAG` is set to the version of the toolkit you want to use.  This will download and add the toolkit to CMake. The toolkit can then be linked against using

.. code:: cmake

   # Link against the IPC Toolkit
   target_link_libraries(${PROJECT_NAME} PUBLIC ipc::toolkit)

where :cmake:`PROJECT_NAME` is the name of your library/binary.

.. tip::
   If your :cmake:`IPC_TOOLKIT_GIT_TAG` is a tag (e.g. ``v1.1.1``), then you can use the :cmake:`FetchContent_Declare` argument :cmake:`GIT_SHALLOW TRUE` to download only a single commit. Otherwise, you should use the default :cmake:`GIT_SHALLOW FALSE`.

Dependencies
------------

**All required dependencies are downloaded through CMake** depending on the build options.

The following libraries are used in this project:

* `Eigen <https://eigen.tuxfamily.org/>`__: linear algebra
* `libigl <https://github.com/libigl/libigl>`__: basic geometry functions and predicates
* `TBB <https://github.com/wjakob/tbb>`__: parallelization
* `Tight-Inclusion <https://github.com/Continuous-Collision-Detection/Tight-Inclusion>`__: correct (conservative) CCD
* `spdlog <https://github.com/gabime/spdlog>`__: logging information

Optional
--------

* `robin-map <https://github.com/Tessil/robin-map>`__: faster hash set/map than :cpp:`std::unordered_set`/:cpp:`std::unordered_map`

  * Enable by using the CMake option :cmake:`IPC_TOOLKIT_WITH_ROBIN_MAP`
  * Enabled by default

* `Abseil <https://abseil.io/>`__: hashing utilities

  * Enable by using the CMake option :cmake:`IPC_TOOLKIT_WITH_ABSEIL`
  * Enabled by default

* `rational-cpp <https://github.io/zfergus/rational-cpp>`__: rational arithmetic used for exact intersection checks

    * Enable by using the CMake option :cmake:`IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION`
    * Requires `GMP <https://gmplib.org/>`__ to be installed at a system level

* `Etienne Vouga's Collision Detection Library <https://github.com/evouga/collisiondetection>`__: inexact CCD

  * Included for comparison with the original IPC library
  * Enable by disabling the CMake option :cmake:`IPC_TOOLKIT_WITH_CORRECT_CCD`
  * Replaces the default Tight-Inclusion CCD

Usage
-----

The main functionality is provided in the ``ipc.hpp`` header. Use the prefix directory ``ipc`` to include all header files (e.g. :cpp:`#include <ipc/ipc.hpp>`).

Unit Tests
----------

We provide unit tests for ensuring the correctness of our algorithmic pieces. To enable the unit tests use the CMake option :cmake:`IPC_TOOLKIT_BUILD_UNIT_TESTS`.

.. _dependencies-1:

Dependencies
^^^^^^^^^^^^

The following are downloaded when unit tests are enabled(:cmake:`IPC_TOOLKIT_BUILD_TESTS`)

* `Catch2 <https://github.com/catchorg/Catch2.git>`__: testing framework
* `finite-diff <https://github.com/zfergus/finite-diff>`__: finite-difference comparisons
* `Nlohman's JSON library <https://github.com/nlohmann/json>`__: loading test data from JSON files