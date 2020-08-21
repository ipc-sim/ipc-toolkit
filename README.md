# IPC Toolkit
**A set of reusable functions to integrate IPC into an existing simulation.**

[![Build status](https://github.com/ipc-sim/ipc-toolkit/workflows/Build/badge.svg?event=push)](https://github.com/ipc-sim/ipc-toolkit/actions?query=workflow%3ABuild+branch%3Amaster+event%3Apush)
[![License](https://img.shields.io/github/license/ipc-sim/ipc-toolkit.svg?color=blue)](https://github.com/ipc-sim/ipc-toolkit/blob/master/LICENSE)

## Integrating IPC Toolkit into your project

#### 1. Add it to CMake

The easiest way to add the toolkit to an existing CMake project is to download
it thorugh CMake. If you do not already have `DownloadProject.cmake` and `DownloadProject.CMakeLists.cmake.in`, you can copy those in our `cmake` folder.
With those added to your project all that is need is to create a file similar to our
`cmake/IPCToolkitDownloadExternal.cmake`. Next, you need to create a file similar to our `cmake/IPCToolkitDependencies.cmake` with the following:

```CMake
# ipc-toolkit
if(NOT TARGET IPCToolkit)
    ${PROJECT_NAME}_download_ipc_toolkit()
    add_subdirectory(${${PROJECT_NAME}_EXTERNAL}/ipc-toolkit EXCLUDE_FROM_ALL)
endif()
```

Lastly, link against the library with

```CMake
target_link_libraries(${PROJECT_NAME} PUBLIC IPCToolkit)
```

#### 2. Use the toolkit

The main functionality is provided in the `ipc.hpp` header. 

## Dependencies

All dependancies are downloaded through CMake depending on the build options.
The following libraries are used in this project:

* [Etienne Vouga's Collision Detection Library](https://github.com/evouga/collisiondetection.git) for continuous collision detection between triangle meshes in 3D
* [libigl](https://github.com/libigl/libigl) as both a source of Eigen and for basic geometry functions (e.g. computing barycentric coordinates)
* [TBB](https://github.com/wjakob/tbb) for limited parallelization
* [Catch2](https://github.com/catchorg/Catch2.git) for testing (see [Unit Tests](#unit_tests))
* [finite-diff](https://github.com/zfergus/finite-diff) for finite difference comparisons in the unit tests

## <a name="unit_tests"></a>Unit Tests

We provide unit tests for ensuring the correctness of our algorithmic pieces.
To enable the unit tests use the flag `-DIPC_TOOLKIT_BUILD_UNIT_TESTS=ON` with
CMake.
