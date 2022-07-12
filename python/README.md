<!-- # IPC Toolkit – Python Bindings -->

[![Python](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/python.yml/badge.svg)](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/python.yml)
[![Docs](https://github.com/ipc-sim/ipc-toolkit/actions/workflows/docs.yml/badge.svg)](https://ipc-sim.github.io/ipc-toolkit/)
[![License](https://img.shields.io/github/license/ipc-sim/ipc-toolkit.svg?color=blue)](https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE)

**A set of reusable functions to integrate IPC into an existing simulation.**

We provide Python bindings for functions in the toolkit using [pybind11](https://github.com/pybind/pybind11).

### Build and Install

Currently, the bindings must be built from scratch. The easiest way to do this is to use the `setup.py` script which uses `setuptools`. For example:
```sh
python setup.py install
```
will build the library and python bindings and then install them on your system.

You can test the install was successful by doing `python -c "import ipctk"`.

### Examples

We provide a Jupyter notebook (`python/example.ipynb`) with some simple examples.

## Contributing

This project is open to contributors! Contributions can come in the form of feature requests, bug fixes, documentation, tutorials and the like. We highly recommend filing an Issue first before submitting a Pull Request.

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
    journal = {ACM Trans. Graph. (SIGGRAPH)},
    year = {2020},
    volume = {39},
    number = {4},
    articleno = {49}
}
```

## License

MIT License © 2020, the IPC-Sim organization (See <a href="https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE"><code>LICENSE.txt</code></a> for details)
