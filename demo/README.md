# Demo Simulator

A small rigid/affine body dynamics simulator built from IPC Toolkit components
(`demo/`, target `ipc_toolkit_demo`, CMake option
`IPC_TOOLKIT_BUILD_DEMO`). It exists to test the components, prototype
ideas, and serve as living proof that the toolkit is easy to integrate — it is
**not** part of the library.

## The library/demo boundary

The core library (`ipc::toolkit`) provides modular, dependency-light math:

- Energy classes (`ipc/dynamics/{rigid,affine}/`): inertial, body-force, and
  orthogonality (ABD) terms, plus the barrier potential — stateless math with
  explicit arguments (`energy(bodies, x)` / `gradient` / `hessian`).
- Time integration (`ipc/dynamics/time_integration/`): BDF (state history +
  `dt` → predicted state); it knows nothing about solvers or stepping loops.
- Contact, friction, adhesion potentials and CCD (unchanged).

This demo (`ipc::toolkit::demo`) owns everything application-side, in two
source pairs:

- `demo/src/simulator.hpp` — the `Simulator`: it *is* a
  [polysolve](https://github.com/polyfem/polysolve)
  `nonlinear::Problem` (energy/gradient/Hessian summed directly from the
  library's energy classes, the toolkit's CCD routed into polysolve's line
  search via `max_step_size`, and the joint-constraint change of variables
  x = Uz [Chen et al. 2022, §4.1] applied around the interface), plus the
  time-stepping loop, adaptive barrier stiffness, pose history, and a stored
  `polysolve::nonlinear::Solver`.
- `demo/src/kinematics.hpp` — how the simulation state maps to the world:
  the rigid (6-DOF) and affine (12-DOF) parameterizations today, with their
  world maps, chain rules, and the mode-appropriate collision candidates and
  CCD. Future discretizations (FEM nodal DOFs, SPH particles) slot in as
  additional `Kinematics` implementations.

The core library never depends on the demo, on polysolve, or on any solver.

## Differences from the original Rigid IPC solver

The solve matches Rigid IPC's [`newton_solver.cpp`](https://github.com/ipc-sim/rigid-ipc/blob/main/src/solvers/newton_solver.cpp)
[Ferguson et al. 2021] in everything that matters — PSD-projected Newton with
an `Eigen::SimplicialLDLT` direction solve, a CCD-filtered backtracking line
search that accepts any energy decrease, and the world-space
max-vertex-velocity stopping criterion (checked from the second iteration,
via the `step_norm()` override) with a gradient-norm early exit. The one
algorithmic difference is the **regularization schedule**: Rigid IPC adapted
a diagonal regularization coefficient continuously inside one Newton loop
(doubled on a failed solve, halved on success, persisting across iterations
and time steps); polysolve instead escalates between discrete strategies —
`ProjectedNewton` → `RegularizedProjectedNewton` (whose internal weight ramps
1e-8 → 1e8) — resetting each solve. Rigid IPC's gradient-descent failsafe is
also dropped: at barrier stiffness scales a −∇f direction is ~1e14 long, so
its CCD-capped line search never yields a usable step.

The remaining non-default polysolve settings in `Settings::solver_params`
(`Backtracking` line search, no plain-Newton stage, `residual_tolerance`
effectively disabled, explicit `SimplicialLDLT`) exist to reproduce the
original behavior — see the comments in `simulator.hpp`.

## Example

```cpp
#include <simulator.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>

using namespace ipc;

// Assemble the scene from rest meshes (V, E, F) and initial poses
std::vector<rigid::Pose> poses = { rigid::Pose(...) };
auto bodies = rigid::RigidBodies::build_from_meshes(
    { V }, { E }, { F }, /*densities=*/{ 1000.0 }, poses);
bodies->planes.emplace_back(Eigen::Vector3d(0, 1, 0), 0); // ground

// Sum the energy terms into an incremental potential and minimize it with
// polysolve's Newton solver, stepping in time with the library's BDF
// integrator and filtering the line search with the library's CCD.
demo::Simulator sim(bodies, poses, /*dt=*/0.01);
sim.run(/*t_end=*/1.0);

// Read back the trajectory (or write it to glTF with ipc/io/write_gltf.hpp
// when the library is built with IPC_TOOLKIT_WITH_GLTF=ON)
std::vector<rigid::Pose> final_poses = sim.rigid_poses();
```

The same scene runs in affine mode (`Settings::body_dynamics = AFFINE`) with
optional joint constraints (`affine::JointConstraints`). The polysolve solver
is configured through `Settings::solver_params` /
`Settings::linear_solver_params` (raw polysolve json; unset entries get the
defaults of polysolve's `nonlinear-solver-spec.json`). polysolve logs through
its own logger, silenced by default (it reports non-gradient stops — e.g.,
every velocity-criterion stop — at error level and per-iteration timings at
debug); set `Settings::solver_log_level` to re-enable it.

## Python

Building with `IPC_TOOLKIT_BUILD_PYTHON=ON` (and the demo enabled) exposes the
simulator as the `ipctk.demo` submodule of the core module:

```python
import ipctk

settings = ipctk.demo.Simulator.Settings()
settings.solver_params = {**settings.solver_params, "max_iterations": 50}
sim = ipctk.demo.Simulator(bodies, initial_poses, dt, settings)
sim.run(1.0)
```

See `python/examples/rigid.py` and `python/examples/joints.py`.

## Tests

The demo's tests live in the unified `ipc_toolkit_tests` binary
(`tests/src/tests/demo/`) and are only compiled when the demo is enabled
(`IPC_TOOLKIT_BUILD_DEMO=ON`). They cover the stepping loop, the joint
constraints, both dynamics models, and finite-difference checks of the
gradient/Hessian through the polysolve `Problem` interface.
