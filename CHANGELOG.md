# Changelog

All notable changes to this project will be documented in this file.

<!--
### YYYY-MM-DD ([XXXXXXX](https://github.com/ipc-sim/ipc-toolkit/commit/XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX))
#### Added
#### Removed
#### Changed
#### Fixed
-->

### 2020-10-07 ([5582582](https://github.com/ipc-sim/ipc-toolkit/commit/5582582bc2f54464bfcee4ba0ec2b7e6975f596f))
#### Added
* `FrictionConstraint` structures to store friction information (i.e., tangent basis, normal force magnitude, closest points, and coefficient of friction)
* Unit test that compares the original IPC code's friction components with the toolkit's.

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
