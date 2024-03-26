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

```
pip install ipctk
```

If you wish to install the current development code, you can compile the library from scratch. Either clone the [repo](https://github.com/ipc-sim/ipc-toolkit) manually or use `git+` with `pip`:

```
pip install git+https://github.com/ipc-sim/ipc-toolkit
```

## Build

To manually build the python binding you can either use the `setup.py` script or use `cmake` directly. The easiest way is to use the `setup.py` script which uses `setuptools`. To do this, use the following command from the root of the repository:

```sh
pip install .
```

This will build the library and install them on your system.

You can test that the installation was successful by doing
```sh
python -c "import ipctk"
```

### CMake Build

Alternatively, you can use `cmake` directly. To do this, use the following commands from the root of the repository:

```sh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DIPC_TOOLKIT_BUILD_PYTHON=ON ..
make -j4
```

## Help/Documentation

* A tutorial on how to use the toolkit can be found [here](https://ipctk.xyz/tutorial/getting_started.html).
* We provide a Jupyter notebook (`python/example.ipynb`) with some simple examples.
* A function reference can be found [here](https://ipctk.xyz/python.html).
