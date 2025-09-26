# [Siggraph 2025] Geometric Contact Potential

<a href="https://github.com/geometryprocessing/GCP-toolkit/actions/workflows/continuous.yml"><img src="https://github.com/geometryprocessing/GCP-toolkit/actions/workflows/continuous.yml/badge.svg"></a>
<a href="https://github.com/geometryprocessing/GCP-toolkit/blob/main/LICENSE"><img src="https://img.shields.io/github/license/geometryprocessing/GCP-toolkit.svg?color=blue"></a>

![Teaser](./docs/teaser.png)

## Description

GCP toolkit is a set of reusable functions to integrate [Geometric Contact Potential](https://huangzizhou.github.io/research/smooth-contact.html) (GCP) into a simulation, built on top of the [IPC toolkit](https://github.com/ipc-sim/ipc-toolkit).

## Features

[Geometric Contact Potential](https://huangzizhou.github.io/research/smooth-contact.html) introduces a new type of barrier-based collision model, with the following features different from IPC:

* **Large $\hat{d}$ allowed**: In IPC, if $\hat{d}$ is larger than any edge length, spurious contact forces would push neighboring vertices apart, which does not happen to GCP.
* **Convergence under refinement**: The GCP potential is derived from a continuous formulation, and the discretized potential converges to the continuous version under mesh refinement.
* **No spurious forces**: GCP elliminates the spurious forces of IPC in various cases (check the paper for details) by considering the local normal and tangent directions.

## Limitations

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

This project is licensed under the MIT License.

You are free to use, modify, and distribute this code in your projects, even commercial ones, as long as you include the original copyright and license notice. A copy of the full license text can be found in the <a href="https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE"><code>LICENSE</code></a> file.

If you use this code in a product you distribute to others, you are required to **include a copy of the original copyright and license notice**. This is typically done in the product's documentation, an "About" or "Third-Party Licenses" section, or in a clear open-source software statement.
