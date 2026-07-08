# Demo Simulator

A small rigid/affine body dynamics simulator built from IPC Toolkit components
(`demo/`, target `ipc_toolkit_demo`, CMake option
`IPC_TOOLKIT_BUILD_DEMO`). It exists to test the components, prototype
ideas, and serve as living proof that the toolkit is easy to integrate — it is
**not** part of the library.

## The library/demo boundary

The core library (`ipc::toolkit`) provides modular, dependency-light math:

- Energy classes (`ipc/dynamics/{rigid,affine}/`): inertial, body-force,
  orthogonality (ABD), and augmented-Lagrangian terms, plus the barrier and
  friction potentials — stateless math with explicit arguments
  (`energy(bodies, x)` / `gradient` / `hessian`). Every term except inertia is
  a **pure, unscaled physical potential** (energy units); the demo applies the
  time-integration scaling when it sums them (see below).
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
  time-stepping loop, adaptive barrier stiffness, the friction/kinematic outer
  loop, pose history, and a stored `polysolve::nonlinear::Solver`.
- `demo/src/kinematics.hpp` — how the simulation state maps to the world:
  the rigid and affine parameterizations (in 2D and 3D), with their world
  maps, chain rules, and the mode-appropriate collision candidates and CCD.
  Future discretizations (FEM nodal DOFs, SPH particles) slot in as additional
  `Kinematics` implementations.

The core library never depends on the demo, on polysolve, or on any solver.

## The incremental potential (objective units)

The step minimizes an incremental potential kept uniformly in **kg·m²**: the
inertial term ½(x − x̃)ᵀM(x − x̃) is unscaled, and **every other term is
multiplied by the integrator's acceleration scaling** `s = (βΔt)²` (= Δt² for
implicit Euler) as the `Simulator` sums them. Because β changes while the BDF
order ramps up, the scaling is applied at evaluation time — never baked into
the term classes or into κ. The library's energy classes are therefore pure
potentials, and the adaptive barrier stiffness balances ∇E against `s·κ·∇B`
(see `initialize_barrier_stiffness`). Barrier stiffness κ is stored only in the
`BarrierPotential`, so the friction term automatically sees the current normal
forces.

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

## Friction

Set `Settings::friction_coefficient > 0` to enable lagged friction
[Ferguson et al. 2021]. Unlike the original's displacement formulation, the
friction potential is evaluated on **velocities** recovered from the time
integrator: `v = (V(x) − Σᵢ wᵢ V(xᵗ⁻ⁱ)) / (βΔt)`, using the same history
combination the BDF scheme uses for its own velocities (so it is correct — not
just consistent — at BDF order ≥ 2). `Settings::static_friction_speed_bound` is
therefore a true speed bound εᵥ (m/s), not a displacement. The tangent bases
and normal-force magnitudes are frozen (lagged) during each inner solve and
rebuilt between solves; `Settings::friction_iterations` controls the lagging:
`0` disables friction, `> 0` caps the number of solves per step (the default
`1` is a single lagged solve), and `< 0` iterates until the momentum balance
‖∇E + κ∇B + ∇D‖ falls below `1e-2 × bbox_diagonal` (or the loop stalls).

> **Accuracy note.** The friction forces are proportional to the *lagged*
> barrier normal forces, so their accuracy follows the accuracy of the contact
> forces at the solved state. Under the default fast velocity-convergence
> criterion those forces are noisy (the barrier's force-vs-gap curve is steep
> near equilibrium). For quantitative friction — matching an analytic
> deceleration μg, for instance — iterate the lagging to convergence
> (`friction_iterations = -1`) **and** tighten `Settings::velocity_conv_tol`
> (e.g. `1e-3`) so the solver cannot stop while the contact forces are
> unbalanced.

## Static and kinematic bodies

Each body has a `type` (`RigidBody::Type::{DYNAMIC, KINEMATIC, STATIC}`,
default `DYNAMIC`):

- **STATIC** bodies never move: all their DOFs are pinned and they contribute
  no inertia, gravity, or external forces (they are still collision-active).
- **KINEMATIC** bodies are driven to a per-step target pose by an augmented
  Lagrangian [Ferguson et al. 2021]. The per-body target policy lives in a
  `KinematicDriver` owned by the `Simulator`: either a scripted absolute pose,
  or (the default) the body's prescribed velocity integrated one step. They
  respond to no forces but do push dynamic bodies through the barrier. A
  driver's optional drive time counts down each step; when it elapses the body
  converts to STATIC.

Set a body's `type` (via `RigidBody::set_type`) before constructing the
`Simulator`; a `KINEMATIC` body gets a default velocity-driven driver, so
`set_type(KINEMATIC)` plus an initial velocity is all a velocity-driven body
needs. To script a body or give it a finite drive time, attach a driver after
construction:

```cpp
(*bodies)[0].set_type(rigid::RigidBody::Type::KINEMATIC);
demo::Simulator sim(bodies, poses, dt, settings);
sim.set_kinematic_driver(0, demo::KinematicDriver::scripted(script));
// or: demo::KinematicDriver::velocity_driven(/*max_time=*/2 * dt)
```

The augmented Lagrangian lives in the library
(`ipc/dynamics/{rigid,affine}::AugmentedLagrangian`) as another energy term
(manual derivatives, no autodiff). The `Simulator` drives it and the friction
lagging in a **single outer loop**: each iteration is one Newton solve followed
by an AL penalty/multiplier update (satisfied channels freeze their DOFs) and a
friction re-lag, until both converge (`Settings::max_outer_iterations` caps it,
with `al_initial_penalty` / `al_max_penalty` / `al_satisfied_progress` /
`al_stall_progress` controlling the schedule).

Individual DOFs of a *dynamic* body can also be fixed by passing an
`is_dof_fixed` mask (`[position | rotation]`, world axes) to
`build_from_meshes`. Joints may not reference a STATIC, KINEMATIC, or
fixed-DOF body (they would be over-constrained); use `add_fixed_point` /
`add_fixed_body` to attach a body to the world instead.

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

// Optionally make a body static or kinematic before constructing the sim
(*bodies)[0].set_type(rigid::RigidBody::Type::STATIC);

// Sum the energy terms into an incremental potential and minimize it with
// polysolve's Newton solver, stepping in time with the library's BDF
// integrator and filtering the line search with the library's CCD.
demo::Simulator::Settings settings;
settings.friction_coefficient = 0.3; // optional lagged friction
demo::Simulator sim(bodies, poses, /*dt=*/0.01, settings);
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
settings.friction_coefficient = 0.3
settings.solver_params = {**settings.solver_params, "max_iterations": 50}
sim = ipctk.demo.Simulator(bodies, initial_poses, dt, settings)
sim.run(1.0)
```

Body classification is exposed on `ipctk.RigidBody` (`type`, `is_dof_fixed`,
`convert_to_static()`); kinematic driving is `ipctk.demo.KinematicDriver`
(`velocity_driven` / `scripted`) attached via
`Simulator.set_kinematic_driver(body, driver)`. See `python/examples/rigid.py`
and `python/examples/joints.py`.

## Tests

The demo's tests live in the unified `ipc_toolkit_tests` binary
(`tests/src/tests/demo/`) and are only compiled when the demo is enabled
(`IPC_TOOLKIT_BUILD_DEMO=ON`). They cover the stepping loop, joint constraints,
both dynamics models, friction (finite-difference derivatives plus analytic
stick/slip and deceleration), static/kinematic bodies (both modes), 2D
simulation, codimensional contact, and finite-difference checks of the
gradient/Hessian through the polysolve `Problem` interface. The
augmented-Lagrangian terms are finite-difference-tested in
`tests/src/tests/dynamics/`.
