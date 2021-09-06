# IPC Toolkit
**A set of reusable functions to integrate IPC into an existing simulation.**

[![Build](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml/badge.svg)](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml)
[![Python](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/python.yml/badge.svg)](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/python.yml)
[![CodeFactor](https://www.codefactor.io/repository/github/ipc-sim/ipc-toolkit/badge)](https://www.codefactor.io/repository/github/ipc-sim/ipc-toolkit)
[![License](https://img.shields.io/github/license/ipc-sim/ipc-toolkit.svg?color=blue)](https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE)

For a complete list of changes, please see [`CHANGELOG.md`](https://github.com/ipc-sim/ipc-toolkit/blob/main/CHANGELOG.md), and for a definitive reference for these functions, please see the [IPC source code](https://github.com/ipc-sim/IPC).

## Integrating IPC Toolkit into your project

#### 1. Add it to CMake

The easiest way to add the toolkit to an existing CMake project is to download it through CMake.
CMake provides functionality for doing this called [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html) (requires CMake ≥ 3.14).
We use this same process to download all external dependencies.
For example,

```CMake
include(FetchContent)
FetchContent_Declare(
    ipc-toolkit
    GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
    GIT_TAG ${IPC_TOOLKIT_GIT_TAG}
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(ipc-toolkit)
```

where `IPC_TOOLKIT_GIT_TAG` is set to the version of the toolkit you want to use. This will download and add the toolkit to CMake. The toolkit can then be linked against using

```CMake
# Link against the IPC Toolkit
target_link_libraries(${PROJECT_NAME} PUBLIC ipc::toolkit)
```

where `PROJECT_NAME` is the name of your project.

#### 2. Using the toolkit

The main functionality is provided in the `ipc.hpp` header. Use the prefix directory `ipc` to include all header files (e.g. `#include <ipc/ipc.hpp>`).

## Dependencies

**All dependancies are downloaded through CMake** depending on the build options.
The following libraries are used in this project:

* [Eigen](https://eigen.tuxfamily.org/): linear algebra
* [libigl](https://github.com/libigl/libigl): basic geometry functions and predicates
* [TBB](https://github.com/wjakob/tbb): parallelization
* [Tight Inclusion CCD](https://github.com/Continuous-Collision-Detection/Tight-Inclusion): correct (conservative) continuous collision detection between triangle meshes in 3D

### Optional

* [Etienne Vouga's Collision Detection Library](https://github.com/evouga/collisiondetection): continuous collision detection between triangle meshes in 3D
    * Enable by using the CMake flag `-DIPC_TOOLKIT_WITH_CORRECT_CCD=OFF`
* [spdlog](https://github.com/gabime/spdlog): logging information (enabled by default)
    * Disable logging completely using the CMake flag `-DIPC_TOOLKIT_WITH_LOGGER=OFF`
* [fmt](https://github.com/fmtlib/fmt): string formatting
    * This is either provided through spdlog or downloaded directly if logging is disabled (`-DIPC_TOOLKIT_WITH_LOGGER=OFF`)
* [Catch2](https://github.com/catchorg/Catch2.git): testing (see [Unit Tests](#unit_tests))
* [finite-diff](https://github.com/zfergus/finite-diff): finite difference comparisons
    * Only used by the unit tests (if they are enabled)


## <a name="unit_tests"></a>Unit Tests

We provide unit tests for ensuring the correctness of our algorithmic pieces.
To enable the unit tests use the flag `-DIPC_TOOLKIT_BUILD_UNIT_TESTS=ON` with
CMake.

## Contributing

This project is open for contributors! Contibutions can come in the form of feature requests, bug fixes, documentation, tutorials and the like. We highly recommend to file an Issue first before submitting a Pull Request.

Simply fork this repository and make a Pull Request! We'd definitely appreciate:

* Implementation of new features
* Bug Reports
* Documentation
* Testing

## License

MIT License © 2020, the IPC-Sim organization (See [`LICENSE.txt`](https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE) for details)
