# Changelog

All notable changes to this project will be documented in this file.

<!--
### YYYY-MM-DD ([XXXXXXX](https://github.com/ipc-sim/ipc-toolkit/commit/XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX))
#### Added
#### Removed
#### Changed
#### Fixed
-->

### 2021-04-29 ([293d0ad](https://github.com/ipc-sim/ipc-toolkit/commit/293d0ad992c01df561e25c286043c9ae9b901ff0))
#### Added
* Added TBB parallel loops to the main function (`compute_potential`, `compute_friction_potential`, `compute_collision_free_stepsize`, etc.)
* Added function `addVerticesFromEdges` that adds the vertices connected to edges in parallel and avoids duplicates

#### Changed
* Changed the HashGrid to use `ArrayMax3` over `VectorX3` to simplify the code

#### Fixed
* Fixed some parameters that were not by reference

### 2021-04-21 ([c8a6d5](https://github.com/ipc-sim/ipc-toolkit/commit/c8a6d56823793e7be5e89238c3793e25bc45ffa0))
#### Added
* Added the SpatialHash from the original IPC code base with some modification to get all candidates in parallel
    * Benchmark results indicate this SpatialHash is faster than the HashGrid with multithreading
    * TODO: Improve HashGrid or fully integrate SpatialHash into ipc.hpp

### 2021-02-11 ([9c7493](https://github.com/ipc-sim/ipc-toolkit/commit/9c74938fefa691db6b79c73489c8c661638019c6))
#### Changed
* Switched to the correct (conservative) CCD of [[Wang and Ferguson et al. 2020]](https://continuous-collision-detection.github.io/)
    * Can select Etienne Vouga's CCD in the CMake (see README.md)

### 2021-02-01 ([b510253](https://github.com/ipc-sim/ipc-toolkit/commit/b51025310223b487e7c39858265d8d5c3e8b1e8a))
#### Added
* Added minimum seperation distance (thickness) to distance constraints
    * Based on [Codimensional Incremental Potential Contact [Li et al. 2020]](https://arxiv.org/abs/2012.04457)

### 2021-02-01 ([a395175](https://github.com/ipc-sim/ipc-toolkit/commit/a3951750ca5f167ab1d546ae1dadd87d0a9e2497))
#### Added
* Added 2D friction model based on the 3D formulation.
    * TODO: Test this further

### 2021-01-12 ([deee6d0](https://github.com/ipc-sim/ipc-toolkit/commit/deee6d0f9802910c5565f800492f9a995e65cf7e))
#### Added
* Added and optional parameter `F2E` to `construct_constraint_set()`. This is similar to `F` (which maps faces to vertices), but maps faces to edges. This is optional, but recommended for better performance. If not provided a simple linear search will be done per face edge!
    * TODO: Add a function to compute this mapping.

### 2021-01-09 ([deee6d0](https://github.com/ipc-sim/ipc-toolkit/commit/deee6d0f9802910c5565f800492f9a995e65cf7e))
#### Changed
* Replaced VectorXd and MatrixXd with static size versions for local gradient and hessians

### 2020-11-20 ([93143ad](https://github.com/ipc-sim/ipc-toolkit/commit/93143ad9b31030cde7324a83354268021e1cb9da))
#### Changed
* Removed TBB parallelization form the hash grid because we get better performance without it.
   * TODO: Improve parallelization in the hash grid or switch to the original IPC spatial hash

### 2020-11-06 ([4553509](https://github.com/ipc-sim/ipc-toolkit/commit/4553509fe6a4e6b78c041018cd6db3fdf23b4730))
#### Fixed
* Fixed multiplicity for point-triangle distance computation to avoid duplicate point-point and point-edge pairs.

### 2020-10-22 ([51f4903](https://github.com/ipc-sim/ipc-toolkit/commit/51f49030dbeec15a6a7544826f5531811a779402))
#### Fixed
* Projection of the hessian to PSD. This was completely broken as the projected matrix was never used.

### 2020-10-22 ([9be6c0f](https://github.com/ipc-sim/ipc-toolkit/commit/9be6c0f7e2534e426e3f09f4c547406d50d5cf9c))
#### Fixed
* Mollification of EE constraints that have a distance type of PP or PE
* If there is no mollification needed then the PP and PE constraints are stored with multiplicity
* Set the parallel EE friction constraint threshold to eps_x like in IPC
    * This avoid needing the mollification for the normal force and these forces are small anyways

### 2020-10-10 ([cb8b53f](https://github.com/ipc-sim/ipc-toolkit/commit/cb8b53fb098598ba5e8c95d4bdb4730e8df9382e))
#### Fixed
* Assertions in `compute_collision_free_stepsize`

### 2020-10-10 ([4a5f84f](https://github.com/ipc-sim/ipc-toolkit/commit/4a5f84f1177bdae1a265dc15a84603bbc389936d))
#### Fixed
* Point-triangle distance type by replacing it with the one used in the original IPC code

### 2020-10-10 ([1d51a61](https://github.com/ipc-sim/ipc-toolkit/commit/1d51a61d60bb25e08c9937285ff9e44459a2223f))
#### Added
* Boolean parameter in `compute_friction_potential_hessian` that controls if the hessian is projected to PSD

### 2020-10-09 ([b737fb0](https://github.com/ipc-sim/ipc-toolkit/commit/b737fb0e708eac5a7775766f162a5d2067db2fa4))
#### Added
* Parameter for vertex group IDs to exclude some collisions (e.g., self collisions)

### 2020-10-08 ([6ee60ae](https://github.com/ipc-sim/ipc-toolkit/commit/6ee60aeaef6d7f88013ee2ee3d544e7403282527))
#### Added
* Second version of `update_barrier_stiffness()` that takes an already computed minimum distance and world bounding box diagonal

### 2020-10-08 ([cc3947d](https://github.com/ipc-sim/ipc-toolkit/commit/cc3947d48bc069488f6a773424e30fe67eb4b5f1))
#### Added
* Second version of `initial_barrier_stiffness()` that takes an already computed barrier gradient
* Assertions on `initial_barrier_stiffness()` input
    * `average_mass > 0 && min_barrier_stiffness_scale > 0`

#### Changed
* Fixed typo in `initial_barrier_stiffness()` name (was `intial_barrier_stiffness()`)

### 2020-10-07 ([5582582](https://github.com/ipc-sim/ipc-toolkit/commit/5582582bc2f54464bfcee4ba0ec2b7e6975f596f))
#### Added
* `FrictionConstraint` structures to store friction information (i.e., tangent basis, normal force magnitude, closest points, and coefficient of friction)
* Unit test that compares the original IPC code's friction components with the toolkit's

#### Changed
* `compute_friction_bases()` is now `construct_friction_constraint_set()`
    * It now takes the coefficient of friction (`mu`)
    * It now puts all information inside of the `FrictionConstraints` (`friction_constraint_set`)

### 2020-10-06 ([b48ba0e](https://github.com/ipc-sim/ipc-toolkit/commit/b48ba0ec9d60754e7670e28fd1987b0c78cd809f))

#### Changed
* During `construct_constraint_set()` the constraints are added based on distance type
    * Duplicate vertex-vertex and edge-vertex constraints are handled by a multiplicity multiplier
    * Edge-edge constraints are always line-line distances
    * Point-triangle constraints are always point-plane distances

### 2020-10-05 ([9a4576b](https://github.com/ipc-sim/ipc-toolkit/commit/9a4576b209302c79296593ac213ed8ce85510f3b))

#### Fixed
* Fixed a bug in the point-triangle closest points and tangent basis computed in `compute_friction_bases()`
* Fixed a bug in `edge_edge_tangent_basis()` used to compute the tangent basis for friction


### 2020-09-19 ([31a37e0](https://github.com/ipc-sim/ipc-toolkit/commit/31a37e04abc9ecec325e00be97fd42b89c895b45))

#### Added
* spdlog for logging information

### 2020-09-19 ([acb7664](https://github.com/ipc-sim/ipc-toolkit/commit/acb7664792982685f6de28468ba126f5e531834f))

#### Changed
* Headers are now include with the prefix `ipc/`
    * E.g., `#include <ipc.hpp>` → `#include <ipc/ipc.hpp>`


### 2020-09-04 ([7dd2ab7](https://github.com/ipc-sim/ipc-toolkit/commit/7dd2ab7a255ffd23ccdfe5aee08bca6a142f75a7))

#### Added
* Collision constraint to store distance constraint pairs
    * `EdgeEdgeConstraint` stores the edge-edge mollifier threshold (`eps_x`)

#### Changed
* Input parameter `dhat_squared` is now `dhat` (i.e., non-squared value)
* Input parameter `epsv_times_h_squared` is now `epsv_times_h` (i.e., non-squared value)
* `Constraints` replaced `Candidates`
* `construct_constraint_set()` now takes the rest vertex position (`V_rest`)
* `compute_barrier_potential*()` no longer take the rest vertex position
