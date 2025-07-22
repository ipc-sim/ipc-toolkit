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

## Build

Instruction for building and including the IPC Toolkit in your CMake project can be found on the website [here](https://ipctk.xyz/build.html).

### Dependencies

The IPC Toolkit depends on a handful of third-party libraries, which are used to provide various functionality.

**All required dependencies are downloaded through CMake** depending on the build options, and are built automatically when you build the IPC Toolkit. You do not need to install them separately.

A full list of dependencies can be found on the [dependencies page](https://ipctk.xyz/dependencies.html).

## Python Bindings

We provide Python bindings for functions in the toolkit using [pybind11](https://github.com/pybind/pybind11).

For more information see the [Python documentation](https://ipctk.xyz/python.html).

## Usage

See the [tutorials](https://ipctk.xyz/tutorials/getting_started.html) for a quick introduction to the toolkit, or the [documentation](https://ipctk.xyz/cpp-api/potentials.html) for a full reference.

## Contributing

This project is open to contributors! Contributions can come in the form of feature requests, bug fixes, documentation, tutorials, and the like. We highly recommend filing an Issue first before submitting a Pull Request.

Simply fork this repository and make a Pull Request! We would appreciate:

* Implementation of new features
* Bug Reports
* Documentation
* Testing

## Citation

IPC Toolkit is created and maintained by academics: citations let us know our work is having impact! Please cite the IPC Toolkit or otherwise give a shout-out if and when it contributes to published works.

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
