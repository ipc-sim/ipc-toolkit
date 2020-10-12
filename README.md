# IPC Toolkit
**A set of reusable functions to integrate IPC into an existing simulation.**

[![Build status](https://github.com/ipc-sim/ipc-toolkit/workflows/Build/badge.svg?event=push)](https://github.com/ipc-sim/ipc-toolkit/actions?query=workflow%3ABuild+branch%3Amaster+event%3Apush)
[![License](https://img.shields.io/github/license/ipc-sim/ipc-toolkit.svg?color=blue)](https://github.com/ipc-sim/ipc-toolkit/blob/master/LICENSE)

:warning: This toolkit is in an early stage of development. If you have any problems or find any bugs please post an issue. For a complete list of changes, please see [`CHANGELOG.md`](https://github.com/ipc-sim/ipc-toolkit/blob/master/CHANGELOG.md). Meanwhile, for a definitive reference for these functions, please see the [IPC source code](https://github.com/ipc-sim/IPC).

## Integrating IPC Toolkit into your project

#### 1. Add it to CMake

The easiest way to add the toolkit to an existing CMake project is to download
it through CMake. If you do not already have `DownloadProject.cmake` and `DownloadProject.CMakeLists.cmake.in`, you can copy those in our `cmake` folder.
With those added to your project, you can create a file a similar to `cmake/IPCToolkitDownloadExternal.cmake` with

```CMake
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
  set(PROJECT_NAME_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
  set(PROJECT_NAME_EXTRA_OPTIONS "")
endif()

function(project_name_download_project name)
  download_project(
    PROJ         ${name}
    SOURCE_DIR   ${PROJECT_NAME_EXTERNAL}/${name}
    DOWNLOAD_DIR ${PROJECT_NAME_EXTERNAL}/.cache/${name}
    QUIET
    ${PROJECT_NAME_EXTRA_OPTIONS}
    ${ARGN}
  )
endfunction()

################################################################################

# Download the IPC Toolkit
function(project_name_download_ipc_toolkit)
  project_name_download_project(ipc-toolkit
    GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
    GIT_TAG        ${IPC_TOOLKIT_GIT_TAG}
  )
endfunction()
```

where `PROJECT_NAME` is the name of your project and `IPC_TOOLKIT_GIT_TAG` is set to the version of the toolkit you want. Then add the toolkit to CMake with

```CMake
# Add the IPC Toolkit to CMake
if(NOT TARGET IPCToolkit)
    project_name_download_ipc_toolkit()
    add_subdirectory(${PROJECT_NAME_EXTERNAL}/ipc-toolkit EXCLUDE_FROM_ALL)
endif()
```

where `PROJECT_NAME_EXTERNAL` is the directory you want the toolkit to be downloaded to. The last step is to link against the library with

```CMake
# Link against the IPC Toolkit
target_link_libraries(${PROJECT_NAME} PUBLIC IPCToolkit)
```

#### 2. Using the toolkit

The main functionality is provided in the `ipc.hpp` header. Use the prefix directory `ipc` to include all header files (e.g. `#include <ipc/ipc.hpp>`).

## Dependencies

All dependancies are downloaded through CMake depending on the build options.
The following libraries are used in this project:

* [Etienne Vouga's Collision Detection Library](https://github.com/evouga/collisiondetection.git) for continuous collision detection between triangle meshes in 3D
* [libigl](https://github.com/libigl/libigl) as both a source of Eigen and for basic geometry functions (e.g. computing barycentric coordinates)
* [TBB](https://github.com/wjakob/tbb) for limited parallelization
* [spdlog](https://github.com/gabime/spdlog) for logging information
    * Disable logging completely using the CMake flag `-IPC_TOOLKIT_WITH_LOGGER=OFF`
* [Catch2](https://github.com/catchorg/Catch2.git) for testing (see [Unit Tests](#unit_tests))
* [finite-diff](https://github.com/zfergus/finite-diff) for finite difference comparisons in the unit tests

## <a name="unit_tests"></a>Unit Tests

We provide unit tests for ensuring the correctness of our algorithmic pieces.
To enable the unit tests use the flag `-DIPC_TOOLKIT_BUILD_UNIT_TESTS=ON` with
CMake.
