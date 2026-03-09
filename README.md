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

* IPC barrier function and its derivatives, and adaptive barrier stiffness algorithm
* Broad- and narrow-phase continuous collision detection (CCD) of linear and nonlinear trajectories
* Distance computation and derivatives between edges in 2D and triangles in 3D
* Distance barrier potential and its derivatives
* Smooth and lagged dissipative friction potential and its derivatives

### Limitations

This is not a full simulation library. As such, it does not include any physics or solvers. For a full simulation implementation, we recommend [PolyFEM](https://polyfem.github.io/) (a finite element library) or [Rigid IPC](https://github.com/ipc-sim/rigid-ipc) (rigid-body dynamics), both of which utilize the IPC Toolkit.

## Build

Instructions for building and including the IPC Toolkit in your CMake project can be found on the website [here](https://ipctk.xyz/build/c++.html).

### Dependencies

The IPC Toolkit depends on a handful of third-party libraries, which are used to provide various functionality.

**All required dependencies are downloaded through CMake** depending on the build options, and are built automatically when you build the IPC Toolkit. You do not need to install them separately.

A full list of dependencies can be found on the [dependencies page](https://ipctk.xyz/build/python.html).

## Python Bindings

We provide Python bindings for functions in the toolkit using [pybind11](https://github.com/pybind/pybind11).

For more information, see the [Python documentation](https://ipctk.xyz/python.html).

## Quick Start

**C++** (include via CMake; see [build instructions](https://ipctk.xyz/build/c++.html)):

```c++
#include <ipc/ipc.hpp>
#include <Eigen/Core>

// Two parallel triangles separated by a small gap (just inside dhat)
const double dhat = 1e-3;
const double gap = 0.5 * dhat; // separation between the two triangles
Eigen::MatrixXd vertices(6, 3);
vertices <<  0.0, 0.0, 0.0,  // triangle 1
             1.0, 0.0, 0.0,
             0.5, 1.0, 0.0,
             0.0, 0.0, gap,  // triangle 2, shifted in z by `gap`
             1.0, 0.0, gap,
             0.5, 1.0, gap;
Eigen::MatrixXi edges(6, 2), faces(2, 3);
edges << 0, 1, 1, 2, 2, 0, 3, 4, 4, 5, 5, 3;
faces << 0, 1, 2, 3, 4, 5;

// Build the collision mesh
ipc::CollisionMesh collision_mesh(vertices, edges, faces);

// Detect collisions
// The two triangles are within dhat, so collisions are active
ipc::NormalCollisions collisions;
collisions.build(collision_mesh, vertices, dhat);

// Compute the barrier potential and its derivatives
const ipc::BarrierPotential B(dhat, /*stiffness=*/1e3);
double energy = B(collisions, collision_mesh, vertices);
Eigen::VectorXd gradient = B.gradient(collisions, collision_mesh, vertices);
Eigen::SparseMatrix<double> hessian = B.hessian(collisions, collision_mesh, vertices);

// Compute the collision-free step size for CCD
Eigen::MatrixXd next_vertices = vertices;
next_vertices.bottomRows(3).col(2).array() -= gap * 2; // move toward triangle 1
double max_step = ipc::compute_collision_free_stepsize(
    collision_mesh, vertices, next_vertices);
```

**Python** (install via `pip install ipctk`):

```python
import ipctk
import numpy as np

# Two parallel triangles separated by a small gap (just inside dhat)
dhat = 1e-3  # barrier activation distance
gap = 0.5 * dhat  # separation between the two triangles
vertices = np.array([
    [0.0, 0.0, 0.0],   # triangle 1
    [1.0, 0.0, 0.0],
    [0.5, 1.0, 0.0],
    [0.0, 0.0, gap],   # triangle 2, shifted in z by `gap`
    [1.0, 0.0, gap],
    [0.5, 1.0, gap],
])
edges = np.array([[0, 1], [1, 2], [2, 0], [3, 4], [4, 5], [5, 3]])
faces = np.array([[0, 1, 2], [3, 4, 5]])

# Build the collision mesh
collision_mesh = ipctk.CollisionMesh(vertices, edges, faces)

# Detect collisions (the two triangles are within dhat, so collisions are active)
collisions = ipctk.NormalCollisions()
collisions.build(collision_mesh, vertices, dhat)

# Compute the barrier potential and its derivatives
B = ipctk.BarrierPotential(dhat, stiffness=1e3)
energy = B(collisions, collision_mesh, vertices)
gradient = B.gradient(collisions, collision_mesh, vertices)
hessian = B.hessian(collisions, collision_mesh, vertices)

# Compute the collision-free step size for CCD
next_vertices = vertices.copy()
next_vertices[3:, 2] -= gap * 2  # move the second triangle toward the first
max_step = ipctk.compute_collision_free_stepsize(
    collision_mesh, vertices, next_vertices)
```

For more details, see the [tutorials](https://ipctk.xyz/tutorials/getting_started.html) and [API documentation](https://ipctk.xyz/cpp-api/potentials.html).

## Contributing

This project is open to contributors! Contributions can come in the form of feature requests, bug fixes, documentation, tutorials, and the like. We highly recommend filing an Issue first before submitting a Pull Request.

Simply fork this repository and make a Pull Request! We would appreciate:

* Implementation of new features
* Bug Reports
* Documentation
* Testing

## Citation

The IPC Toolkit is created and maintained by academics: citations let us know our work is having an impact! Please cite the IPC Toolkit or otherwise give a shout-out if and when it contributes to published works.

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
