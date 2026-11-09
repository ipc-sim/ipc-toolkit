# Python Bindings

![PyPI](https://img.shields.io/pypi/v/ipctk?color=brightgreen&label=PyPI&logo=python&logoColor=white)
![PyPI - Downloads](https://img.shields.io/pypi/dm/ipctk?label=PyPI%20Downloads&logo=python&logoColor=white)
[![Python](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/python.yml/badge.svg)](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/python.yml)
[![Docs](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/docs.yml/badge.svg)](https://ipctk.xyz/)
[![codecov](https://codecov.io/github/ipc-sim/ipc-toolkit/graph/badge.svg?token=9BR6GPKRY8)](https://codecov.io/github/ipc-sim/ipc-toolkit)
[![License](https://img.shields.io/github/license/ipc-sim/ipc-toolkit.svg?color=blue)](https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE)

We provide Python bindings for functions in the toolkit using [pybind11](https://github.com/pybind/pybind11).

## Installation

To install the latest release, you can use `pip`:

```sh
pip install ipctk
```

If you wish to install the current development code, you can compile the library from scratch. Either clone the [repo](https://github.com/ipc-sim/ipc-toolkit) manually or use `git+` with `pip`:

```sh
pip install git+https://github.com/ipc-sim/ipc-toolkit
```

## Build

To manually build the python binding you can either use `pip` or use `cmake` directly. The easiest way is `pip`, which builds the bindings with [scikit-build-core](https://scikit-build-core.readthedocs.io). To do this, use the following command from the root of the repository:

```sh
pip install .
```

This will build the library and install them on your system.

You can test that the installation was successful by doing
```sh
python -c "import ipctk"
```

The build lives in `build/<wheel tag>/`, so re-running `pip install .` only recompiles what changed.

#### CMake options

Pass CMake options with `-C cmake.define.<NAME>=<VALUE>`, or through the `CMAKE_ARGS` environment variable:

```sh
pip install . -C cmake.define.IPC_TOOLKIT_WITH_CUDA=ON
CMAKE_ARGS="-DIPC_TOOLKIT_WITH_CUDA=ON" pip install .
```

To build in debug mode, set `SKBUILD_CMAKE_BUILD_TYPE=Debug` (or pass `-C cmake.build-type=Debug`).

#### SIMD

SIMD optimizations are enabled by default. To disable them (e.g., for compatibility with older or cross-compiled targets), set `IPCTK_WITH_SIMD=0` before installing:

```sh
IPCTK_WITH_SIMD=0 pip install .
```

Accepted falsy values: `0`, `off`, `false`, `no` (case-insensitive). Any other value keeps SIMD enabled.

:::{note}
Pre-built binary wheels from PyPI have SIMD support baked in at wheel-build time. The `IPCTK_WITH_SIMD` variable only applies when building from source. To force a source build from PyPI, use:
```sh
IPCTK_WITH_SIMD=0 pip install --no-binary ipctk ipctk
```
:::

### CMake Build

Alternatively, you can use `cmake` directly. To do this, use the following commands from the root of the repository:

```sh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DIPC_TOOLKIT_BUILD_PYTHON=ON ..
make -j$(nproc)
```

## Help/Documentation

* Tutorials on how to use the toolkit can be found [here](https://ipctk.xyz/tutorials/getting_started.html).
* We provide a Jupyter notebook (`python/examples/example.ipynb`) with some simple examples.
* A function reference can be found [here](https://ipctk.xyz/).
