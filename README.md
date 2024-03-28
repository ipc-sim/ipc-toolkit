<p align="center">
<a href="https://ipctk.xyz"><img alt="IPC Toolkit" src="docs/source/_static/logo.png" width="80%"></a>
</p>

<p align="center">
<a href="https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml"><img src="https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml/badge.svg"></a>
<a href="https://github.com/ipc-sim/ipc-toolkit/actions/workflows/python.yml"><img src="https://github.com/ipc-sim/ipc-toolkit/actions/workflows/python.yml/badge.svg"></a>
<a href="https://ipctk.xyz"><img src="https://github.com/ipc-sim/ipc-toolkit/actions/workflows/docs.yml/badge.svg"></a>
<a href="https://codecov.io/github/ipc-sim/ipc-toolkit"><img src="https://codecov.io/github/ipc-sim/ipc-toolkit/graph/badge.svg?token=9BR6GPKRY8"/></a>
<a href="https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE"><img src="https://img.shields.io/github/license/ipc-sim/ipc-toolkit.svg?color=blue"></a>
</p>

## Description

IPC Toolkit is a set of reusable functions to integrate Incremental Potential Contact (IPC) into a simulation.

### Features

* IPC barrier function and its derivatives and adaptive barrier stiffness algorithm
* Broad- and narrow-phase continuous collision detection (CCD) of linear and nonlinear trajectories
* Distance computation and derivatives between edges in 2D and triangles in 3D
* Distance barrier potential and its derivatives
* Smooth and lagged dissipative friction potential and its derivatives

### Limitations

This is not a full simulation library. As such it does not include any physics or solvers. For a full simulation implementation, we recommend [PolyFEM](https://polyfem.github.io/) (a finite element library) or [Rigid IPC](https://github.com/ipc-sim/rigid-ipc) (rigid-body dynamics) both of which utilize the IPC Toolkit.

<!--- BEGIN C++ README 1 --->

## Build

The easiest way to add the toolkit to an existing CMake project is to download it through CMake.
CMake provides functionality for doing this called [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html) (requires CMake ≥ 3.14).
We use a very similar process to download all external dependencies (using [CPM](https://github.com/cpm-cmake/CPM.cmake)).

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

<!--- BEGIN C++ README 2 --->

### Dependencies

**All required dependencies are downloaded through CMake** depending on the build options.

The following libraries are used in this project:

* [Eigen](https://eigen.tuxfamily.org/): linear algebra
* [libigl](https://github.com/libigl/libigl): basic geometry functions and predicates
* [oneTBB](https://github.com/oneapi-src/oneTBB): parallelism
* [Tight-Inclusion](https://github.com/Continuous-Collision-Detection/Tight-Inclusion): provably conservative CCD of [Wang and Ferguson et al. 2021]
* [SimpleBVH](https://github.com/ipc-sim/SimpleBVH): a simple bounding volume hierarchy data structure
* [Scalable-CCD](https://github.com/Continuous-Collision-Detection/Scalable-CCD): scalable (GPU) CCD of [Belgrod et al. 2023]
* [spdlog](https://github.com/gabime/spdlog): logging information

#### Optional

The following dependencies are optionally used based on CMake options:

* [robin-map](https://github.com/Tessil/robin-map): faster hash set/map than `std::unordered_set`/`std::unordered_map`
    * Enable by using the CMake option `IPC_TOOLKIT_WITH_ROBIN_MAP`
    * Enabled by default
* [Abseil](https://abseil.io/): hashing utilities
    * Enable by using the CMake option `IPC_TOOLKIT_WITH_ABSEIL`
    * Enabled by default
* [filib](https://github.com/zfergus/filib): interval arithmetic for nonlinear trajectories/CCD
    * Enable by using the CMake option `IPC_TOOLKIT_WITH_FILIB`
    * Enabled by default
* [rational-cpp](https://github.io/zfergus/rational-cpp): rational arithmetic used for exact intersection checks
    * Enable by using the CMake option `IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION`
    * Requires [GMP](https://gmplib.org/) to be installed at a system level
* [Etienne Vouga's Collision Detection Library](https://github.com/evouga/collisiondetection): inexact CCD
    * Included for comparison with the original IPC library
    * Enable by using the CMake option `IPC_TOOLKIT_WITH_INEXACT_CCD`
    * Replaces the default Tight-Inclusion CCD

## Usage

See the [tutorial](https://ipctk.xyz/tutorial/getting_started.html) for a quick introduction to the toolkit, or the [documentation](https://ipctk.xyz/cpp.html) for a full reference.

## Unit Tests

We provide unit tests to ensure the correctness of our algorithmic pieces.
To enable the unit tests use the CMake option `IPC_TOOLKIT_BUILD_TESTS`.

### Dependencies

The following are downloaded when unit tests are enabled:

* [Catch2](https://github.com/catchorg/Catch2.git): testing framework
* [finite-diff](https://github.com/zfergus/finite-diff): finite-difference comparisons
* [Nlohman's JSON library](https://github.com/nlohmann/json): loading test data from JSON files

<!--- END C++ README --->

## Python Bindings

We provide Python bindings for functions in the toolkit using [pybind11](https://github.com/pybind/pybind11).

For more information see the [Python documentation](https://ipctk.xyz/python.html).

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
  url = {https://github.com/ipc-sim/ipc-toolkit},
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

MIT License © 2020, the IPC-Sim organization (See <a href="https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE"><code>LICENSE</code></a> for details).
