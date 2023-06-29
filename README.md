# IPC Toolkit

[![Build](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml/badge.svg)](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml)
[![Python](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/python.yml/badge.svg)](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/python.yml)
[![Docs](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/docs.yml/badge.svg)](https://ipc-sim.github.io/ipc-toolkit/)
[![License](https://img.shields.io/github/license/ipc-sim/ipc-toolkit.svg?color=blue)](https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE)

## Description

IPC Toolkit is a set of reusable functions to integrate Incremental Potential Contact (IPC) into a simulation.

### Features

* IPC barrier function and its derivatives and adaptive barrier stiffness algorithm
* Broad-phase and narrow-phase continuous collision detection (CCD)
* Distance computation and derivatives between edges in 2D and triangles in 3D
* Distance barrier potential and its derivatives
* Smooth and lagged dissipative friction potential and its derivatives

### Limitations

This is not a full simulation library. As such it does not include any physics or solvers. For a full simulation implementation, we recommend [PolyFEM](https://polyfem.github.io/) (a finite element library) or [Rigid IPC](https://github.com/ipc-sim/rigid-ipc) (rigid-body dynamics) both of which utilize the IPC Toolkit.

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

#### Optional

* [robin-map](https://github.com/Tessil/robin-map): faster hash set/map than `std::unordered_set`/`std::unordered_map`
    * Enable by using the CMake option `IPC_TOOLKIT_WITH_ROBIN_MAP`
    * Enabled by default
* [Abseil](https://abseil.io/): hashing utilities
    * Enable by using the CMake option `IPC_TOOLKIT_WITH_ABSEIL`
    * Enabled by default
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

## Python Bindings

We provide Python bindings for functions in the toolkit using [pybind11](https://github.com/pybind/pybind11).

For more information see the [Python documentation](https://ipc-sim.github.io/ipc-toolkit/python/).

## Contributing

This project is open to contributors! Contributions can come in the form of feature requests, bug fixes, documentation, tutorials, and the like. We highly recommend filing an Issue first before submitting a Pull Request.

Simply fork this repository and make a Pull Request! We would appreciate:

* Implementation of new features
* Bug Reports
* Documentation
* Testing

## Citation

If you use the IPC Toolkit in your project, please consider citing our work:

```bibtex
@software{ipc_toolkit,
  author = {Zachary Ferguson and others},
  title = {{IPC Toolkit}},
  url = {https://ipc-sim.github.io/ipc-toolkit/},
  year = {2020},
}
```

Additionally, you can cite the original IPC paper:

```bibtex
@article{Li2020IPC,
    author = {Minchen Li and Zachary Ferguson and Teseo Schneider and Timothy Langlois and
        Denis Zorin and Daniele Panozzo and Chenfanfu Jiang and Danny M. Kaufman},
    title = {Incremental Potential Contact: Intersection- and Inversion-free Large Deformation Dynamics},
    journal = {{ACM} Trans. Graph. (SIGGRAPH)},
    year = {2020},
    volume = {39},
    number = {4},
    articleno = {49}
}
```

## License

MIT License © 2020, the IPC-Sim organization (See <a href="https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE"><code>LICENSE</code></a> for details)
