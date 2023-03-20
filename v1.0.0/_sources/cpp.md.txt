# C++

[![Build](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml/badge.svg)](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml)
[![Docs](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/docs.yml/badge.svg)](https://ipc-sim.github.io/ipc-toolkit/)
[![License](https://img.shields.io/github/license/ipc-sim/ipc-toolkit.svg?color=blue)](https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE)

## Build

The easiest way to add the toolkit to an existing CMake project is to download it through CMake.
CMake provides functionality for doing this called [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html) (requires CMake ≥ 3.14).
We use this same process to download all external dependencies.

For example,

```cmake
include(FetchContent)
FetchContent_Declare(
    ipc_toolkit
    GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
    GIT_TAG ${IPC_TOOLKIT_GIT_TAG}
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(ipc_toolkit)
```

where `IPC_TOOLKIT_GIT_TAG` is set to the version of the toolkit you want to use. This will download and add the toolkit to CMake. The toolkit can then be linked against using

```cmake
# Link against the IPC Toolkit
target_link_libraries(${PROJECT_NAME} PUBLIC ipc::toolkit)
```

where `PROJECT_NAME` is the name of your library/binary.

### Dependencies

**All required dependencies are downloaded through CMake** depending on the build options.

The following libraries are used in this project:

* [Eigen](https://eigen.tuxfamily.org/): linear algebra
* [libigl](https://github.com/libigl/libigl): basic geometry functions and predicates
* [TBB](https://github.com/wjakob/tbb): parallelization
* [Tight-Inclusion](https://github.com/Continuous-Collision-Detection/Tight-Inclusion): correct (conservative) CCD
* [spdlog](https://github.com/gabime/spdlog): logging information
* [robin-map](https://github.com/Tessil/robin-map): faster hash set/map than `std::unordered_set`/`std::unordered_map`
* [Abseil](https://abseil.io/): hashing utilities

#### Optional

* [GMP](https://gmplib.org/): rational arithmetic used for exact intersection checks
    * Enable by using the CMake option `IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION`
    * GMP must be installed at a system level
* [Etienne Vouga's Collision Detection Library](https://github.com/evouga/collisiondetection): inexact CCD
    * Included for comparison with the original IPC library
    * Enable by disabling the CMake option `IPC_TOOLKIT_WITH_CORRECT_CCD`
    * Replaces the default Tight-Inclusion CCD

## Usage

The main functionality is provided in the `ipc.hpp` header. Use the prefix directory `ipc` to include all header files (e.g. `#include <ipc/ipc.hpp>`).

## Unit Tests

We provide unit tests for ensuring the correctness of our algorithmic pieces.
To enable the unit tests use the CMake option `IPC_TOOLKIT_BUILD_UNIT_TESTS`.

### Dependencies

The following are downloaded when unit tests are enabled (`IPC_TOOLKIT_BUILD_TESTS`)

* [Catch2](https://github.com/catchorg/Catch2.git): testing framework
* [finite-diff](https://github.com/zfergus/finite-diff): finite-difference comparisons
* [Nlohman's JSON library](https://github.com/nlohmann/json): loading test data from JSON files