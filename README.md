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

This is not a full simulation library. As such it does not include any physics or solvers. A full simulation implementation in [PolyFEM](https://polyfem.github.io) is available [here](https://github.com/geometryprocessing/geometric-contact-potential).

## Usage

[Geometric Contact Potential](https://huangzizhou.github.io/research/smooth-contact.html) differs from IPC only in the potential formulation. It reuses existing utilities in IPC toolkit to collect collision candidates, perform narrow phase CCD, etc. Below shows a table of GCP classes with their IPC counterparts.

| GCP | IPC | Description |
| ----- | ------ | ----- |
| SmoothCollision | Collision | Base class of various collision pair type |
| SmoothCollisionTemplate | EdgeEdgeCollision, FaceVertexCollision, etc. | Specific collision pair type |
| SmoothCollisionsBuilder | CollisionsBuilder | To construct collision pairs from collision candidates |
| SmoothCollisions | Collisions | Collection of all collision pairs |
| SmoothContactPotential | BarrierPotential | Collision potential |