# Changelog

All notable changes to this project will be documented in this file.

## v1.1.1 (Aug 18, 2023)

* Logo by @zfergus in [#52](https://github.com/ipc-sim/ipc-toolkit/pull/52)
* Fix vertex-vertex `==` and `<` functions to be order independent
    * This allows vertex-vertex constraints merge correctly
* Update Tight Inclusion CCD

## v1.1.0 (Jul 25, 2023)

Large refactoring to make the code more object-oriented rather than passing objects to functions. Other changes include the friction potential now being a function of velocity, bug fixes, and a new tutorial.

### Details

* Large Refactor in [#25](https://github.com/ipc-sim/ipc-toolkit/pull/25)
    * `construct_collision_candidates(..., candidates)` -> `candidates.build(...)`
    * `is_step_collision_free(candidates, ...)` -> `candidates.is_step_collision_free(...)`
    * `compute_collision_free_stepsize(candidates, ...)` -> `candidates.compute_collision_free_stepsize(...)`
    * `compute_barrier_potential*(constraints, ...)` -> `constraints.compute_potential*(...)`
    * `compute_shape_derivative(constraints, ...)` -> `constraints.compute_shape_derivative(...)`
    * `compute_minimum_distance(constraints, ...)` -> `constraints.compute_minimum_distance(...)`
    * `construct_friction_constraint_set(..., friction_constraints)` -> `friction_constraints.build(...)`
    * `compute_friction_*(..., friction_constraints, ...)` -> `friction_constraints.compute_*(...)`
    * Generic `CollisionStencil` parent class to `Candidates`, `CollisionConstraints`, and `FrictionConstraints`.
    * Renamed `Constraints` to `CollisionConstraints`
    * Replaced single letter variable names `V`, `E`, `F` with `vertices`/`positions`, `edges`, `faces`
    * Renamed `*_index` -> `*_id`
    * Replaced `inflation_radius = min_distance / 1.99` with `inflation_radius = min_distance / 2` and use rounding mode to conservativly inflate AABBs
    * `CollisionConstraints::use_convergent_formulation` and `are_shape_derivatives_enabled` must now be accessed through getter and setter functions
    * Friction potentials are now functions of velocity. Previously `V0` and `V1` were passed and `U = V1-V0`. This limited the integration scheme to implicit Euler. Upstream this means you need to multiply the potential by `1/(dv/dx)` to get the correct friction force.
        * Change input $\epsilon_vh$ to $\epsilon_v$ in [#37](https://github.com/ipc-sim/ipc-toolkit/pull/37) to reflect the fact that friction is defined in terms of velocity instead of displacement now.
* Changed default `project_hessian_to_psd` to `false` in [#30](https://github.com/ipc-sim/ipc-toolkit/pull/30)
* Update website with a tutorial ([#31](https://github.com/ipc-sim/ipc-toolkit/pull/31)) and version dropdown list ([#34](https://github.com/ipc-sim/ipc-toolkit/pull/34))
* Switch from templates to using `Eigen::Ref` in [#28](https://github.com/ipc-sim/ipc-toolkit/pull/28)
* Speed up the CCD by limiting the maximum minimum distance to `1e-4` in [#43](https://github.com/ipc-sim/ipc-toolkit/pull/43)
* Fix the bug pointed out in [#41](https://github.com/ipc-sim/ipc-toolkit/pull/41) in [#42](https://github.com/ipc-sim/ipc-toolkit/pull/42). Namely, to get units of distance in the barrier we should divide the original function by $\hat{d}\cdot(\hat{d} + 2d_{\min})^2$ when using distance squared. Before it was being divided by $2d_{\min} \hat{d} + \hat{d}^2$.
* Fix build for IPC_TOOLKIT_WITH_CORRECT_CCD=OFF in [#44](https://github.com/ipc-sim/ipc-toolkit/pull/44)
* Switched from FetchContent to CPM in [#48](https://github.com/ipc-sim/ipc-toolkit/pull/48). This provides better caching between builds. Additionally, made robin-map and Abseil optional dependencies.
* Add the CFL-Inspired Culling of CCD as described in Section 3 of the Technical Supplement to IPC in [#50](https://github.com/ipc-sim/ipc-toolkit/pull/50)

## v1.0.0 (Feb 21, 2023)

This is the first official release. 🚀

This is a stable release of the toolkit prior to refactoring the code and making updates to the API.

### Details

* Added a minimum distance optional parameter to all CCD functions (`const double min_distance = 0.0`) in [#22](https://github.com/ipc-sim/ipc-toolkit/pull/22). This is placed as the first optional argument which can break calling code if optional parameters were previously used.
* Added `CollisionMesh` in [#7](https://github.com/ipc-sim/ipc-toolkit/pull/7) to wrap up face and edges into a single data structure.
    * Removes Support for ignoring internal vertices. Instead, users should use the CollisionMesh to map from the full mesh to the surface mesh.
    * This also includes a `to_full_dof` function that can map the reduced gradient/hessian to the full mesh's DOF.

## Pre-v1.0.0

### 2021-10-05 ([9e2cc2a](https://github.com/ipc-sim/ipc-toolkit/commit/574f7577daa5e0b51bf5baf20998994b8371216e))
#### Added
* Added implicits source folder to organize point-plane collisions

### 2021-09-05 ([9e2cc2a](https://github.com/ipc-sim/ipc-toolkit/commit/9e22cc2a5f7e7ca048a579f2c94d2241782ecf17))
#### Added
* Added support for point vs. (static) analytical plane contact

### 2021-08-21 ([acf2a80](https://github.com/ipc-sim/ipc-toolkit/commit/acf2a80544ebe27dc5e440602a3a89243e575e8a))
#### Changed
* Changed CMake target name to `ipc::toolkit`

### 2021-07-26 ([1479aae](https://github.com/ipc-sim/ipc-toolkit/commit/1479aaea958daaa4e963529493e4169dc7757913))
#### Changed
* Updated the CMake system to use modern `FetchContent` to download externals

### 2021-07-22 ([e24c76d](https://github.com/ipc-sim/ipc-toolkit/commit/e24c76ddc818fb9efc4d522ef72a581a15abf751))
#### Fixed
* Updated CCD strategy when using Tight Inclusion to only perform `no_zero_toi=true` when there is no minimum distance

### 2021-07-17 ([a20f7a2](https://github.com/ipc-sim/ipc-toolkit/commit/a20f7a2dfea5a04c67ef71d0cd523f69391f2f54))
#### Added
* Added `detect_edge_face_collision_candidates_brute_force` for 3D intersection broad-phase
* Added ability to save an obj of collision candidates
* Added tests for has_intersection (all pass after fixes)

#### Fixed
* Fixed possible numerical rounding problems in HashGrid `AABB::are_overlapping`
* Fixed HashGrid's function for getting edge-face intersection candidates

### 2021-07-15 ([7301b42](https://github.com/ipc-sim/ipc-toolkit/commit/7301b422a9b9a90c76d9e7abf2f9127bf6d0dbd6))
#### Fixed
* Use `ignore_codimensional_vertices` in the brute force broad-phase method
* Fixed AABB inflation in brute force and SpatialHash methods

### 2021-07-08 ([86ae4e5](https://github.com/ipc-sim/ipc-toolkit/commit/86ae4e5f87eb2c65585920ad3ca0bbb3b57702f6))
#### Changed
* Replaced vertex group ids with more powerful can_collide function. By
default everything can collide with everything (same as before)
* Reordered parameters in `construct_constraint_set()`,
`is_collision_free()`, and `compute_collision_free_stepsize()`
* `update_barrier_stiffness` now requires the `constraint_set` rather
than building it
* `update_barrier_stiffness` dropped dhat parameter

#### Fixed
* SpatialHash for 2D

#### Removed
* Verison of `initial_barrier_stiffness` that computes the
constraint set and barrier gradient because there are a lot of
parameters to these functions

### 2021-07-05 ([4d16954](https://github.com/ipc-sim/ipc-toolkit/commit/4d16954012570b3a15346b99b5aedea77266fe86))
#### Changed
* Renamed directory `src/spatial_hash/` → `src/broad_phase/`
* Renamed files `src/ccd/broad_phase.*` → `src/ccd/aabb.*`

### 2021-07-05 ([b3808e1](https://github.com/ipc-sim/ipc-toolkit/commit/b3808e15bdbaba9a6efd4b731db3070e85bcc4b7))
#### Added
* Select the broad-phase method for CCD and distance constraints
    * Methods: `HASH_GRID`, `SPATIAL_HASH`, `BRUTE_FORCE`
* CCD parameters for Tight Inclusion's tolerance and maximum iterations

#### Changed
* `ignore_codimensional_vertices` to `false` by default
* CMake option `TIGHT_INCLUSION_WITH_NO_ZERO_TOI=ON` as default

### 2021-06-18 ([aa59aeb](https://github.com/ipc-sim/ipc-toolkit/commit/aa59aeb0634af981a8f1cfbb6d2ff2b76a04d610))
#### Changed
* `construct_friction_constraint_set` now clears the given `friction_constraint_set`

### 2021-05-18 ([245b13b](https://github.com/ipc-sim/ipc-toolkit/commit/245b13bcc5e99ed52850ae865aaa0ad4e71a43a8))
#### Changed
* Use TightInclusion degenerate edge-edge for point-point and point-edge CCD

### 2021-05-11 ([5c34dcd](https://github.com/ipc-sim/ipc-toolkit/commit/5c34dcdf226d46ada962204585fa386eb9b67859))
#### Changed
* `char*` exceptions to `std::exceptions`

### 2021-05-06 ([24056cc](https://github.com/ipc-sim/ipc-toolkit/commit/24056ccb2ca0a03bdef8141bc5011c41547f06b5))
#### Changed
* Gave `dhat_epsilon_scale` a default value of `1e-9` in `update_barrier_stiffness`
* :warning: Changed order of parameters to `update_barrier_stiffness`
    * Flipped `bbox_diagonal` and `dhat_epsilon_scale`

### 2021-05-06 ([81d65f3](https://github.com/ipc-sim/ipc-toolkit/commit/81d65f32e479fea32d0acc29c8a7a532fa55518b))
#### Fixed
* Bug in output min distance of `update_barrier_stiffness`

### 2021-05-04 ([59ec167](https://github.com/ipc-sim/ipc-toolkit/commit/59ec167b85eaf56095a2d0333bdd96146d658ebf))
#### Changed
* Moved eigen_ext functions into ipc namespace
* Renamed max size matrices with `Max`
    * `Eigen::VectorX([0-9])` → `ipc::VectorMax$1`
    * `Eigen::MatrixXX([0-9])` → `ipc::VectorMax$1`
    * `Eigen::ArrayMax([0-9])` → `ipc::ArrayMax$1`

### 2021-05-03 ([664d65f](https://github.com/ipc-sim/ipc-toolkit/commit/664d65fd70dbd350b6bfe5f8a311a89ff4fef3bd))
#### Added
* Added utility function to check for edge-edge intersection in 2D and edge-triangle intersection in 3D.
* Optionally: use GMP for exact edge-triangle intersection checks

### 2021-05-03 ([9b4ebfc](https://github.com/ipc-sim/ipc-toolkit/commit/9b4ebfc0f458645cf33eeebf8211607f45ad9cb4))
#### Added
* voxel_size_heuristic.cpp which suggests a good voxel size for the SpatialHash and HashGrid

#### Changed
* Changed HashGrid voxel size to be the average edge length not considering displacement length. This results in better performance, but can result in large memory usage.

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
