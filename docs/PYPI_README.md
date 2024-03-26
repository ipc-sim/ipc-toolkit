# IPC Toolkit

![PyPI](https://img.shields.io/pypi/v/ipctk?color=brightgreen&label=PyPI&logo=python&logoColor=white)
![PyPI - Downloads](https://img.shields.io/pypi/dm/ipctk?label=PyPI%20Downloads&logo=python&logoColor=white)
![GitHub Repo stars](https://img.shields.io/github/stars/ipc-sim/ipc-toolkit?label=Stars&logo=github)
[![codecov](https://codecov.io/github/ipc-sim/ipc-toolkit/graph/badge.svg?token=9BR6GPKRY8)](https://codecov.io/github/ipc-sim/ipc-toolkit)
[![License](https://img.shields.io/github/license/ipc-sim/ipc-toolkit.svg?color=blue&label=License)](https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE)

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

## Installation

To install the latest release, you can use `pip`:

```
pip install ipctk
```

If you wish to install the current development code, you can compile the library from scratch. Either clone the [repo](https://github.com/ipc-sim/ipc-toolkit) manually or use `git+` with `pip`:

```
pip install git+https://github.com/ipc-sim/ipc-toolkit
```

## Help/Documentation

* A tutorial on how to use the toolkit can be found [here](https://ipctk.xyz/tutorial/getting_started.html).
* A function reference can be found [here](https://ipctk.xyz/python.html).

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
    journal = {ACM Trans. Graph. (SIGGRAPH)},
    year = {2020},
    volume = {39},
    number = {4},
    articleno = {49}
}
```

## License

MIT License © 2020, the IPC-Sim organization (See <a href="https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE"><code>LICENSE.txt</code></a> for details)
