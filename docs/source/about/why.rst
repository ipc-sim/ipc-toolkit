Why the IPC Toolkit
===================

The IPC Toolkit is a set of functions for computing contact in a simulation. Instead of adopting a new simulator, you call these functions from your own solver, time integrator, and data structures.

Contact is the hardest, most failure-prone part of any physics simulation. Getting it intersection-free, with correct derivatives, robust continuous collision detection (CCD), and stable friction takes years of research and engineering. Many people already have a solver and a mesh pipeline that they are happy with. What they are missing is contact that they can trust, and that is the gap the toolkit fills.

.. note::
    The IPC Toolkit is **not** a full simulation library. It does not include any physics or solvers. We think of that as the point of the project rather than a limitation: we implement the contact and leave everything else to you.

What You Get
------------

.. list-table::
    :header-rows: 1
    :widths: 25 75

    * - Component
      - What it provides
    * - Barrier potential
      - The IPC log-barrier with its gradient, Hessian (optionally projected to be positive semi-definite), and the adaptive barrier stiffness algorithm. See :doc:`/tutorials/getting_started`.
    * - Continuous collision detection
      - Broad- and narrow-phase CCD for linear and nonlinear trajectories, and a collision-free step size to use in your line search. See :doc:`/tutorials/nonlinear_ccd`.
    * - Distances
      - Edge–edge distances in 2D and triangle-based distances in 3D, with derivatives.
    * - Friction and adhesion
      - Smooth, lagged dissipative friction and adhesion potentials with derivatives. See :doc:`/tutorials/advanced_friction` and :doc:`/tutorials/adhesion`.
    * - Alternative formulations
      - The convergent formulation, the Geometric Contact Potential, and Offset Geometric Contact, alongside the original IPC formulation. See :doc:`/tutorials/convergent`, :doc:`/tutorials/gcp`, and :doc:`/tutorials/ogc`.
    * - Languages
      - Modern C++17 built with CMake, and Python bindings with ``pip install ipctk``.
    * - License
      - MIT, so you can use it in commercial products. See :doc:`license`.

Where the Toolkit Fits
----------------------

Simulation software generally falls into one of two groups. Either you build on top of it, or you embed it inside your own solver.

- **Integrated frameworks** own the whole stack: the scene description, solver, integrator, and I/O. You get a complete simulator, but you have to work within its architecture.
- **Component libraries** solve one problem and hand the result back to you. You keep control of everything else.

Both are useful. However, as far as we know, the projects with IPC-grade contact are full frameworks that you build on, and the projects that you can embed do not handle contact. The IPC Toolkit aims to do both.

.. raw:: html
   :file: ../_static/landscape.html

*Illustrative positioning based on each project's stated scope. This is not a benchmark, and it leaves out many other projects. The horizontal axis asks whether you embed the project in your solver or build on top of it. The vertical axis asks how closely the project's contact handling resembles IPC (a barrier with intersection-free CCD).*

The table below lists the projects in the figure, starting with the ones closest to the IPC Toolkit. Compare the *Embed in your solver?* and *IPC-grade contact?* columns. Projects with IPC-grade contact are almost always frameworks that you build on, and the projects you can embed in your own solver almost never handle contact.

.. list-table::
    :header-rows: 1
    :widths: 16 34 18 18 14

    * - Project
      - What it is
      - Embed in your solver?
      - IPC-grade contact?
      - Actively maintained?
    * - **IPC Toolkit**
      - Modular IPC contact, friction, and CCD components (C++ and Python)
      - Yes, by design
      - Yes
      - Yes
    * - `libuipc <https://github.com/spiriMirror/libuipc>`_
      - Unified GPU IPC framework for rigid, soft, cloth, and thread simulation
      - No, you build on it
      - Yes
      - Yes
    * - `GIPC <https://github.com/KemengHuang/GPU_IPC>`_
      - GPU-accelerated IPC with Gauss–Newton optimization of the barrier energy :cite:`Huang2024GIPC`
      - No, it is its own simulator
      - Yes
      - No, a reference implementation
    * - `IPC <https://github.com/ipc-sim/IPC>`_ and `C-IPC <https://github.com/ipc-sim/Codim-IPC>`_
      - The reference implementations of :cite:t:`Li2020IPC` and :cite:t:`Li2021CIPC`
      - No, it is research code built around one solver
      - Yes
      - No, a reference implementation
    * - `SymX <https://github.com/InteractiveComputerGraphics/SymX>`_
      - C++ symbolic differentiation with JIT code generation and global gradient and Hessian assembly for nonlinear optimization. It is the engine of `STARK <https://github.com/InteractiveComputerGraphics/stark>`_.
      - Partly, you bring your own data but build the simulator in SymX
      - Yes, as IPC energies you write in SymX (see its IPC example and STARK)
      - Yes
    * - `YASPS <https://github.com/txstc55/yasps>`_
      - Symbolic framework that generates CUDA code for extensible, high-performance IPC simulation, mixing free and affine-body parameterizations (`arXiv:2605.23088 <https://arxiv.org/abs/2605.23088>`_)
      - No, you build on it
      - Yes
      - Yes
    * - `STARK <https://github.com/InteractiveComputerGraphics/stark>`_
      - C++ and Python platform for strongly coupled rigid and deformable simulation (FEM, shells, rods, and joints) with frictional contact, built on SymX
      - No, you build on it
      - Yes
      - Yes
    * - `PolyFEM <https://polyfem.github.io/>`_
      - Full finite element simulator
      - No, it is an application
      - Yes, via the IPC Toolkit
      - Yes
    * - `Rigid IPC <https://github.com/ipc-sim/rigid-ipc>`_
      - Intersection-free rigid body dynamics :cite:`Ferguson2021RigidIPC`
      - No, it is an application
      - Yes, via the IPC Toolkit
      - No, a reference implementation
    * - `SfePy <https://github.com/sfepy/sfepy>`_, `FEDOO <https://3mah.github.io/fedoo-docs>`_, `MeshFEM <https://github.com/MeshFEM/MeshFEM>`_
      - Finite element libraries (see the :doc:`gallery`)
      - No, they are frameworks
      - Yes, via the IPC Toolkit
      - Yes
    * - `Tight-Inclusion <https://github.com/Continuous-Collision-Detection/Tight-Inclusion>`_, `Scalable-CCD <https://github.com/Continuous-Collision-Detection/Scalable-CCD>`_
      - Provably conservative CCD, on the CPU and GPU. The IPC Toolkit builds on both.
      - Yes
      - CCD only, no barrier or friction
      - Yes
    * - `FCL <https://github.com/flexible-collision-library/fcl>`_, `Coal <https://github.com/humanoid-path-planner/hpp-fcl>`_
      - Collision and distance queries for robotics and planning
      - Yes
      - Queries only, no contact potentials
      - Yes
    * - `TinyAD <https://github.com/patr-schm/TinyAD>`_, `Bartels <https://github.com/dilevin/Bartels>`_
      - Automatic differentiation and finite element energy collections
      - Yes
      - No
      - Rarely updated
    * - `SimKit <https://github.com/otmanon/simkit>`_
      - Python toolbox for prototyping physics-based animation research
      - Yes
      - No
      - Yes
    * - `Ipopt <https://coin-or.github.io/Ipopt/>`_
      - General nonlinear optimizer using an interior-point method
      - Yes, you write the physics
      - No
      - Yes
    * - `FEBio <https://febio.org/>`_, `Abaqus <https://www.3ds.com/products/simulia/abaqus>`_, `Ansys <https://ansys.synopsys.com/>`_, `LS-DYNA <https://lsdyna.ansys.com/>`_
      - General-purpose finite element packages with their own contact algorithms
      - No, they are applications
      - A different approach
      - Yes
    * - `SOFA <https://www.sofa-framework.org/>`_
      - Interactive simulation framework with a plugin architecture
      - Partly, through plugins
      - A different approach
      - Yes
    * - `Chrono <https://projectchrono.org/>`_, `Drake <https://drake.mit.edu/>`_
      - Multibody dynamics and robotics toolboxes
      - No, they are their own engines
      - A different approach
      - Yes
    * - `Bullet <https://pybullet.org/>`_, `PhysX <https://developer.nvidia.com/physx-sdk>`_, `MuJoCo <https://mujoco.org/>`_
      - Real-time rigid-body and robotics physics engines
      - No, they are their own engines
      - A different approach
      - Yes
    * - `NVIDIA Newton <https://developer.nvidia.com/newton-physics>`_
      - GPU differentiable robotics engine built on Warp
      - No, it is its own engine
      - Not its focus
      - Yes

.. note::
    IPC, C-IPC, Rigid IPC, and GIPC are reference implementations of their papers. They are valuable for reproducing results and studying the algorithms, but the code is not actively maintained, so we would not build a long-lived project on top of them. The IPC Toolkit is maintained by authors of the original IPC paper, and it aims to be the maintained home for the building blocks that these methods share. "Rarely updated" means the project still works but sees infrequent changes.

.. note::
    SymX and YASPS take a different route to a similar goal. You write the energies, including IPC contact, in a symbolic framework, and it generates the derivatives and assembly code for you. The IPC Toolkit instead provides hand-derived, tested contact functions that you call from your own code. The tradeoff is flexibility in how you define energies versus keeping your existing solver and data structures.

.. note::
    "A different approach" is not a criticism. Penalty, augmented Lagrangian, complementarity, and soft-contact models work well for many problems. By construction, though, they do not guarantee an intersection-free state the way a barrier with CCD does, and that guarantee is what we mean by IPC-grade contact.

Why Use a Component Library
---------------------------

Bring Your Own Solver
~~~~~~~~~~~~~~~~~~~~~

Every potential provides an energy, a gradient, and a sparse Hessian (optionally projected to be positive semi-definite). You add these to the objective of whatever Newton, quasi-Newton, or projected solver you already use. There is no framework to adopt and no scene format to learn. See :doc:`/tutorials/simulation` for how the pieces fit into an optimization-based time integrator.

It Solves the Hard Part, and Only the Hard Part
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Writing a Newton loop is a weekend project. Getting intersection-free contact with robust CCD, a collision-free line-search step, and stable friction derivatives is a research program. The IPC Toolkit gives you the research program and leaves the weekend to you.

A Thin Dependency
~~~~~~~~~~~~~~~~~

The core library is centered on Eigen. CMake downloads the dependencies for you, and anything beyond the core (SIMD, faster hashing, interval arithmetic, etc.) is behind an opt-in or opt-out build flag, so you only pay for what you use. See :doc:`dependencies` for the full list and their licenses.

Easy to Wire In
~~~~~~~~~~~~~~~

The toolkit asks very little of your solver. You evaluate the value, gradient, and Hessian, and ask for the maximum feasible step size before your line search. This maps onto almost any optimization-based integrator, so integrating the toolkit usually only requires a short adapter instead of restructuring your code. For example, the demo simulator described in `Where It's Going`_ plugs the toolkit into the nonlinear ``Problem`` interface of `PolySolve <https://github.com/polyfem/polysolve>`_ and passes our collision-free step size to PolySolve's line search.

Maintained and Tested
~~~~~~~~~~~~~~~~~~~~~

The toolkit comes from authors of the original IPC paper :cite:p:`Li2020IPC`, and we continue to maintain it. It has continuous integration, code coverage, published documentation, and Python bindings that match the C++ API. It is meant to be production-friendly and publication-friendly at the same time, so if you use it in your work, please cite it (see :doc:`license`).

Where It's Going
----------------

.. note::
    This section describes the ``dynamics`` branch, which is a big part of the upcoming v2.0.0 release. It is not in the current release yet.

Version 2.0.0 extends the toolkit beyond contact, but we are not turning it into a simulator. Something only belongs in the core if it is dependency-light math and does not make the library an application.

- **Rigid and affine body energies.** Inertia, body forces, the orthogonality potential for rigid bodies, an augmented Lagrangian, and joints are all energy terms that you add as needed. We write each one once in affine coordinates and use a change of variables to apply them to rigid bodies, so rigid and affine bodies share a single formulation.
- **A swappable BDF time integrator.** It takes a state history and a time step and returns the predicted state used by the inertial term. It knows nothing about solvers or stepping loops, so you can replace it with your own integrator.
- **A separate demo simulator.** The ``demo/`` directory combines these components with `PolySolve <https://github.com/polyfem/polysolve>`_, steps them with BDF, and exports the result. It is application code, so the library itself never depends on it, PolySolve, or any other solver.

The "Few Lines of Code" Test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The demo is also a test of our interfaces. If connecting the toolkit's components to an off-the-shelf solver only takes a few lines of code, then the interfaces are doing their job. This makes the demo both an integration test of the components and an example of how to adopt the toolkit.

When to Choose Something Else
-----------------------------

Knowing your lane is part of being useful, so here is when we would point you elsewhere:

- If you need a complete finite element application, `PolyFEM <https://polyfem.github.io/>`_ already ships one, powered by the IPC Toolkit.
- If you are training robot policies at scale, `NVIDIA Newton <https://developer.nvidia.com/newton-physics>`_ is the natural home.
- If you want a turnkey GPU simulator for cloth, soft bodies, and rigid bodies, `libuipc <https://github.com/spiriMirror/libuipc>`_ is built for that.

But if you already have a solver and need contact that you can trust inside it, the IPC Toolkit is for you.

Next Steps
----------

- :doc:`/build/c++` or :doc:`/build/python` to install the toolkit.
- :doc:`/tutorials/getting_started` for a first collision mesh, barrier potential, and CCD query.
- :doc:`/tutorials/simulation` to see how the toolkit slots into a time integrator.
- :doc:`gallery` for projects and papers that already use the toolkit.
