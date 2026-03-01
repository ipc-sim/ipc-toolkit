Release Notes
=============

.. role:: cpp(code)
   :language: c++
.. role:: cmake(code)
   :language: cmake

v1.5.0 (Febuary 5, 2026)
------------------------

Highlights
~~~~~~~~~~

- Implement "Geometric Contact Potential" :cite:t:`Huang2025GCP` by `@Huangzizhou <https://github.com/Huangzizhou>`_ in `#191 <https://github.com/ipc-sim/ipc-toolkit/pull/191>`_
- Implement "Offset Geometric Contact" :cite:t:`Chen2025Offset` by `@zfergus <https://github.com/zfergus>`_ in `#192 <https://github.com/ipc-sim/ipc-toolkit/pull/192>`_
- Implement fully parallel LBVH Broad Phase by `@zfergus <https://github.com/zfergus>`_ in `#188 <https://github.com/ipc-sim/ipc-toolkit/pull/188>`_
- Add geometry utilities for collision normals, signed distances, and dihedral angles with gradients and Hessians.

New Formulations
~~~~~~~~~~~~~~~~

- **Geometric Contact Potential (GCP)** by `@Huangzizhou <https://github.com/Huangzizhou>`_ in `#191 <https://github.com/ipc-sim/ipc-toolkit/pull/191>`_

  - The implementation includes the collision and friction potential with gradients and Hessians.
  - Tutorial available `here <https://ipctk.xyz/tutorials/gcp.html>`_.

- **Offset Geometric Contact (OGC)** by `@zfergus <https://github.com/zfergus>`_ in `#192 <https://github.com/ipc-sim/ipc-toolkit/pull/192>`_

  - OGC is a penetration-free contact model that replaces continuous collision detection (CCD) with a trust-region-based approach.
  - **New Namespace:** Introduced ``ipc::ogc`` to encapsulate OGC-specific functionality.
  - **Trust Region Implementation:**

    - Added ``TrustRegion`` class (``trust_region.hpp``, ``trust_region.cpp``) to manage per-vertex conservative bounds.
    - Implemented ``warm_start_time_step`` to initialize the trust region and handle initial predictions.
    - Implemented ``filter_step`` to scale optimization steps ensuring vertices stay within safe bounds.
    - Implemented ``update`` and ``update_if_needed`` to dynamically resize trust regions based on motion.

  - **Feasible Region Logic:**

    - Added ``feasible_region.hpp`` and ``feasible_region.cpp`` containing geometric predicates (e.g., ``check_vertex_feasible_region``, ``is_edge_edge_feasible``) to verify if primitives are within valid non-penetrating regions.
    - Integrated feasible region checks into the ``NormalCollisions`` class to filter out invalid collision candidates. Enabled via ``set_collision_set_type(NormalCollisions::CollisionSetType::OGC)``.

  - **Step Scaling:** Unlike the original paper's projection method, ``TrustRegion::filter_step`` scales the descent direction $\beta$ to keep vertices on the trust region boundary. This preserves the descent direction, ensuring compatibility with line-search-based solvers.
  - Tutorial available `here <https://ipctk.xyz/tutorials/ogc.html>`_.

Broad Phase
~~~~~~~~~~~

- Parallel CPU LBVH implementation by `@zfergus <https://github.com/zfergus>`_ in `#188 <https://github.com/ipc-sim/ipc-toolkit/pull/188>`_

  - New broad phase collision detection method with memory optimizations on the AABB structure.
  - LBVH Broad Phase:

    - Implemented ``LBVH`` class in ``src/ipc/broad_phase/lbvh.hpp``, providing a Linear Bounding Volume Hierarchy method for broad phase collision detection.
    - Fully parallel build and traversal routines using Intel TBB.
    - SIMD optimizations for AABB overlap tests.
    - Faster than the existing BVH method and all other broad phase methods in IPC Toolkit for large scenes.

      - More than 3x faster to build.
      - Up to 1.5x faster for candidate detection.
      - Detailed performance charts available below.

    - Added ``python/examples/lbvh.py`` to demonstrate usage and visualization.

  - ⚡ Enhancements

    - Performance Optimizations:

      - Added ``DefaultInitAllocator`` to reduce overhead when allocating large arrays of POD types.
      - Replaced ``std::vector<AABB>`` with ``std::vector<AABB, DefaultInitAllocator<AABB>>``.
      - Memory Optimization:

        - Switched ``AABB::min`` and ``AABB::max`` from ``ArrayMax3d`` to ``Eigen::Array3d``.
        - Result: Reduces ``sizeof(AABB)`` from 76 to 64 bytes, allowing AABB to fit into a single cache line (assuming 64-byte lines).

      - Replaced hardware rounding in ``AABB::conservative_inflation`` with software rounding using ``std::nextafter``.

  - 💥 Breaking Changes

    - API & Structural Changes:

      - Moved dimension variable: Removed dim from ``AABB`` and moved it to ``BroadPhase`` to handle 2D/3D data more cleanly. ``dim`` is now set in ``BroadPhase::build``.
      - Handle 2D data by setting AABB's z-components to zero.
      - Constructor Removal: Removed the AABB default constructor and the initialization of ``AABB::vertex_ids``.
      - Type Change: Changed ``to_3D`` in ``ipc/utils/eigen_ext.hpp`` to operate on ``Eigen::Array`` types instead of ``Eigen::Vector``.
      - Deprecated ``BVH`` class in favor of the new ``LBVH`` class for better performance on large scenes.

- Refactor BroadPhase to Build from Vertex Boxes Directly in `#187 <https://github.com/ipc-sim/ipc-toolkit/pull/187>`_
- Rename ``sweep_and_tiniest_queue.cpp`` to ``.cu`` file by `@iiiian <https://github.com/iiiian>`_ in `#203 <https://github.com/ipc-sim/ipc-toolkit/pull/203>`_
- Fix vertex id assignment by `@iiiian <https://github.com/iiiian>`_ in `#204 <https://github.com/ipc-sim/ipc-toolkit/pull/204>`_

  - Fix vertex id assignment logic in ``SweepAndTiniestQueue``.

- Replace optional ``const shared_ptr<BroadPhase>&`` parameters with ``BroadPhase*`` in `#205 <https://github.com/ipc-sim/ipc-toolkit/pull/205>`_

  - Change ``make_default_broad_phase`` to return a ``unique_ptr``

- Pass parameters to Scalable CCD Narrow-Phase by `@antoinebou12 <https://github.com/antoinebou12>`_ in `#208 <https://github.com/ipc-sim/ipc-toolkit/pull/208>`_

  - Replace hardcoded ``TightInclusionCCD::DEFAULT_TOLERANCE`` and ``TightInclusionCCD::DEFAULT_MAX_ITERATION`` with ``narrow_phase.tolerance`` and ``narrow_phase.tolerance`` if passed a ``TightInclusionCCD`` object.
  - If ``narrow_phase_ccd`` is not a ``TightInclusionCCD`` instance, we fall back to default values to preserve backward compatibility.

- Add SIMD support via ``xsimd`` library in `#207 <https://github.com/ipc-sim/ipc-toolkit/pull/207>`_

  - Cross-platform SIMD (Single Instruction, Multiple Data) support using the ``xsimd`` library in the LBVH code.
  - Update the ``IPC_TOOLKIT_WITH_SIMD`` to be enabled by default, and improve CMake logic to detect SIMD capabilities and configure the build accordingly.
  - Disable Eigen's internal vectorization to avoid crashes on Linux.

Math and Geometry Utilities
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Add normal computation and Jacobian methods for collision candidates in `#182 <https://github.com/ipc-sim/ipc-toolkit/pull/182>`_
- Collect normal utilities into ``geometry/normal`` files in `#197 <https://github.com/ipc-sim/ipc-toolkit/pull/197>`_
- Organize geometry and math utilities into ``geometry`` and ``math`` folders in `#198 <https://github.com/ipc-sim/ipc-toolkit/pull/198>`_
- Add signed distance computations for point-plane, point-line, and line-line geometries in `#206 <https://github.com/ipc-sim/ipc-toolkit/pull/206>`_
- Optimized edge_edge_mollifier.cpp with FMA operator in gradient calculations.

Miscellaneous
~~~~~~~~~~~~~

- Updated dependencies:

  - Bump Abseil from ``20250512.1`` to ``20260107.0``
  - Bump Catch2 from ``v3.8.1`` to ``v3.12.0``
  - Bump Eigen from ``3.4.0`` to ``5.0.1``
  - Bump JSON from ``v3.11.2`` to ``v3.12.0``
  - Bump libigl from ``89267b4a80b1904de3f6f2812a2053e5e9332b7e`` to ``v2.6.0``
  - Bump TBB from ``v2022.1.0`` to ``v2022.3.0``
  - Bump Pybind11 from ``v2.13.1`` to ``v3.0.1``
  - Bump robin-map from ``v1.4.0`` to ``v1.4.1``
  - Bump spdlog from ``v1.15.3`` to ``v1.17.0``

- Make most dependencies private

  - Added a figure for default dependencies of the ``ipc::toolkit`` library: https://ipctk.xyz/about/dependencies.html
  - Refactored ``BVH`` class to use unique_ptr to hide ``SimpleBVH``
  - Updated Sweep and Prune and Sweep and Tiniest Queue classes to use a Pimpl pattern for encapsulation.
  - Mark ``tsl::robin_map`` and ``absl::hash`` as private dependencies in CMake. Update the dependency documentation accordingly.
  - Convert ``SpatialHash`` to PImpl idiom to hide ``tsl::robin_map`` and ``absl::hash`` from the public API.

- Add Clang Tidy check

  - Added a ``.clang-tidy`` configuration file with a comprehensive set of checks, exclusions, naming conventions, and options to enforce modern C++ practices and treat all warnings as errors.
  - Introduced a GitHub Actions workflow (``.github/workflows/clang-tidy-check.yml``) to automatically run Clang-Tidy on relevant files for every push and pull request, improving code quality and review automation.
- Simplify CMake CUDA setup by `@iiiian <https://github.com/iiiian>`_ in `#201 <https://github.com/ipc-sim/ipc-toolkit/pull/201>`_
  - Disable CUDA builds by default.
  - Bump CMake minimum version to ``3.24`` to support ``CMAKE_CUDA_ARCHITECTURES="native"``

- Add performance profiling utilities in `#188 <https://github.com/ipc-sim/ipc-toolkit/pull/188>`_

  - Added a (optional) profiler utility (``src/ipc/utils/profiler.hpp``) to measure/record metrics with CSV output support.
  - Added Nlohmann JSON as a dependency for profiler storage.
  - Use ``IPC_TOOLKIT_WITH_PROFILER`` CMake option to enable/disable profiling.

- Refactor to use ``Eigen::ConstRef`` for vertices in `#200 <https://github.com/ipc-sim/ipc-toolkit/pull/200>`_

v1.4.0 (July 22, 2025)
----------------------

Highlights
~~~~~~~~~~

- Add adhesion potentials for sticky interactions, allowing for more realistic simulations of contact between objects.

  - Based on the paper "Augmented Incremental Potential Contact for Sticky Interactions" by :cite:t:`Fang2023AugmentedStickyInteractions`

- Implement new barrier functions from the literature.

  - :cpp:class:`ipc::NormalizedClampedLogBarrier`: unit-less clamped log barrier function.
  - :cpp:class:`ipc::ClampedLogSqBarrier`: clamped log barrier with a quadratic log term from :cite:`Huang2024GIPC`
  - :cpp:class:`ipc::CubicBarrier`: cubic barrier function from :cite:`Ando2024Cubic`

- Add support for separate static and kinetic coefficients in tangential collisions.

What's Changed
~~~~~~~~~~~~~~

- Adhesion in `#85 <https://github.com/ipc-sim/ipc-toolkit/pull/85>`_

  - Implement the adhesion potentials from "Augmented Incremental Potential Contact for Sticky Interactions" :cite:t:`Fang2023AugmentedStickyInteractions`
  - Files moved from `ipc/friction <https://github.com/ipc-sim/ipc-toolkit/tree/v1.3.1/src/ipc/friction>`_ to generic `ipc/tangential <https://github.com/ipc-sim/ipc-toolkit/tree/v1.4.0/src/ipc/tangential>`_ and `ipc/collisions/tangential <https://github.com/ipc-sim/ipc-toolkit/tree/v1.4.0/src/ipc/collisions/tangential>`_ folders.
  - Renamed ``Collision`` classes to :cpp:class:`ipc::NormalCollision`
  - Renamed ``FrictionCollision`` classes to :cpp:class:`ipc::TangentialCollision`
  - Renamed ``DistanceBasedPotential`` to :cpp:class:`ipc::NormalPotential`

    - Added :cpp:func:`ipc::NormalPotential::force_magnitude` and :cpp:func:`ipc::NormalPotential::force_magnitude_gradient` methods

  - Added :cpp:class:`ipc::TangentialPotential` parent class to :cpp:class:`ipc::FrictionPotential` and :cpp:class:`ipc::TangentialAdhesionPotential`
  - `Tutorial <https://ipctk.xyz/tutorials/adhesion.html>`_ added by `@antoinebou12 <https://github.com/antoinebou12>`_ in `#144 <https://github.com/ipc-sim/ipc-toolkit/pull/144>`_
  - Use normal potential in tangential collisions by `@maxpaik16 <https://github.com/maxpaik16>`_ in `#142 <https://github.com/ipc-sim/ipc-toolkit/pull/142>`_

    - Replaces the ``const BarrierPotential&`` parameter with a ``const NormalPotential&``
    - This allows tangential adhesion to be properly set up in simulation.

- Add :cpp:class:`ipc::NormalizedClampedLogBarrier` in `#143 <https://github.com/ipc-sim/ipc-toolkit/pull/143>`_
- Add :cpp:member:`ipc::AdditiveCCD::max_iterations` parameter in `#145 <https://github.com/ipc-sim/ipc-toolkit/pull/145>`_
- Add cubic barrier with elasticity-inclusive dynamic stiffness in `#148 <https://github.com/ipc-sim/ipc-toolkit/pull/148>`_

  - Implement the cubic barrier and elasticity-inclusive dynamic stiffness from :cite:`Ando2024Cubic`
  - Implement :math:`log^2`-barrier from GIPC :cite:`Huang2024GIPC`
  - Add :cpp:func:`ipc::CollisionStencil::compute_coefficients` -- computes the coefficients :math:`\vec{c}` s.t. :math:`d(x) = \Vert \sum c_i x_i\Vert^2`
  - Add :cpp:func:`ipc::CollisionStencil::compute_distance` taking :math:`V`, :math:`E`, :math:`F` directly

- Replace ``BroadPhaseMethod`` enum with ``shared_ptr<BroadPhase>`` in `#152 <https://github.com/ipc-sim/ipc-toolkit/pull/152>`_

  - Python BroadPhase in `#155 <https://github.com/ipc-sim/ipc-toolkit/pull/155>`_

- Remove ``ContinuousCollisionCandidate`` inheritance in `#156 <https://github.com/ipc-sim/ipc-toolkit/pull/156>`_

  - Replace it with the :cpp:class:`ipc::CollisionStencil` class.

- Replace ``const Eigen::Matrix&`` with ``Eigen::ConstRef<Matrix>`` in `#157 <https://github.com/ipc-sim/ipc-toolkit/pull/157>`_

  - This allows for more efficient passing of matrices to functions.

- Separate Static and Kinetic Coefficients of Friction in `#177 <https://github.com/ipc-sim/ipc-toolkit/pull/177>`_

  - Enhancements to tangential collision handling by distinguishing between static and kinetic friction coefficients (``mu_s`` and ``mu_k``) and implementing smooth transitions between them.
  - Updated :cpp:class:`ipc::TangentialCollision` class to replace ``mu`` with separate ``mu_s`` (static friction coefficient) and ``mu_k`` (kinetic friction coefficient). This change applies to the class definition and all related methods.
  - New :cpp:func:`ipc::smooth_mu` functions: Added functions to compute smooth transitions between static and kinetic friction coefficients, including derivatives and related mathematical formulations. These functions improve the modeling of friction forces at varying tangential velocities.
  - Integration with adhesion module: Incorporated smooth friction coefficient calculations into the adhesion module.
  - Tutorial added in `#178 <https://github.com/ipc-sim/ipc-toolkit/pull/178>`_

Python Specific
~~~~~~~~~~~~~~~

- Add faster ``can_collide`` functions to Python in `#135 <https://github.com/ipc-sim/ipc-toolkit/pull/135>`_
- Add ``PyBroadPhase`` and ``PyNarrowPhaseCCD`` classes to wrap the :py:class:`ipctk.BroadPhase` and :py:class:`ipctk.NarrowPhaseCCD` classes, respectively, allowing for custom implementations of these classes in Python.
- Add new constructors to candidate classes (:py:class:`ipctk.EdgeEdgeCandidate`, :py:class:`ipctk.EdgeFaceCandidate`, :py:class:`ipctk.EdgeVertexCandidate`, :py:class:`ipctk.FaceFaceCandidate`, :py:class:`ipctk.FaceVertexCandidate`, :py:class:`ipctk.VertexVertexCandidate`) that accept tuples for easier initialization.
- Include a call to ``std::atexit`` to reset the thread limiter upon program exit to ensure proper cleanup.
- :WARNING: Python :py:class:`ipctk.NarrowPhaseCCD` implementations will not work with multi-threading because of GIL locking.
- Switch from ``py::arg`` to ``_a`` literals in `#168 <https://github.com/ipc-sim/ipc-toolkit/pull/168>`_

Miscellaneous
~~~~~~~~~~~~~

- Add ``/Wall`` and ``/MP`` compiler flags for MSVC in `#134 <https://github.com/ipc-sim/ipc-toolkit/pull/134>`_
- Enable CUDA by default if CUDA is detected by CMake in `#134 <https://github.com/ipc-sim/ipc-toolkit/pull/134>`_
- Devcontainer by `@antoinebou12 <https://github.com/antoinebou12>`_ in `#136 <https://github.com/ipc-sim/ipc-toolkit/pull/136>`_
- Update dependencies for CMake 4.0 in `#158 <https://github.com/ipc-sim/ipc-toolkit/pull/158>`_
- Organize folder structure for generated Xcode project in `#159 <https://github.com/ipc-sim/ipc-toolkit/pull/159>`_

  - Refactor CMake configurations and improve project structure
  - Updated CMake recipes for third-party libraries to set appropriate folder properties for IDE organization.
  - Changed the CPMAddPackage references for scalable-ccd to a more recent commit.
  - Enhanced logging functionality by adding a new ``log_and_throw_error`` function to improve error handling.
  - Improved documentation in the style guide for naming conventions.
  - Move generated ``config.hpp`` file from ``src/ipc`` to ``build/include/src/ipc``
  - Update CUDA builds in ``CMakePresets.json``
  - Remove ``virtual`` tag from :cpp:func:`ipc::BroadPhase::detect_collision_candidates`

- Index typedef in `#161 <https://github.com/ipc-sim/ipc-toolkit/pull/161>`_

  - Add a ``index_t`` typedef for vertex/edge/face ids. Changes default from ``long`` to ``int32_t``.

- Update Tight Inclusion package version to 1.0.6 in `#165 <https://github.com/ipc-sim/ipc-toolkit/pull/165>`_
- Fix cmake for non-top-level inclusion by `@sin3point14 <https://github.com/sin3point14>`_ in `#167 <https://github.com/ipc-sim/ipc-toolkit/pull/167>`_
- Improve and clean up documentation:

  - Fix Error in Docs in `#169 <https://github.com/ipc-sim/ipc-toolkit/pull/169>`_
  - Update docs in `#173 <https://github.com/ipc-sim/ipc-toolkit/pull/173>`_
  - Update docs.yml in `#174 <https://github.com/ipc-sim/ipc-toolkit/pull/174>`_
  - Refactor and clean-up C++ API docs in `#175 <https://github.com/ipc-sim/ipc-toolkit/pull/175>`_


v1.3.1 (Nov 08, 2024)
---------------------

-  Move test data to an external repository in `#113 <https://github.com/ipc-sim/ipc-toolkit/pull/113>`__

   -  Download test data at compile time as needed using CMake

-  Individual convergent formulation flags in `#120 <https://github.com/ipc-sim/ipc-toolkit/pull/120>`__

   -  Replace the ``use_convergent_formulation`` flag in ``Collisions`` with two flags: ``use_area_weighting`` and ``use_improved_max_approximator``.
   -  Move the physical barrier rescaling from collision weights into the ``BarrierPotential`` using a flag ``use_physcial_barrier``.

-  Replace scalar ``tbb::enumerable_thread_specific`` with ``tbb::parallel_reduce`` in `#121 <https://github.com/ipc-sim/ipc-toolkit/pull/121>`__
-  Fix missing ``max_iterations`` and ``tolerance`` variables in ``compute_collision_free_stepsize`` when using ``SWEEP_AND_TINIEST_QUEUE`` by `@antoinebou12 <https://github.com/antoinebou12>`__ in `#123 <https://github.com/ipc-sim/ipc-toolkit/pull/123>`__
-  Add ``cuda.yml`` to test if the code compile with CUDA enabled in `#125 <https://github.com/ipc-sim/ipc-toolkit/pull/125>`__
-  Update filib to allow shared library build in `#122 <https://github.com/ipc-sim/ipc-toolkit/pull/122>`__
-  Fix ``friction_collisions.build`` tutorial (Issue `#126 <https://github.com/ipc-sim/ipc-toolkit/pull/126>`__) in `#127 <https://github.com/ipc-sim/ipc-toolkit/pull/127>`__
-  Sort includes using clang-format in `#129 <https://github.com/ipc-sim/ipc-toolkit/pull/129>`__
-  Add frequently asked questions page to the tutorial documentation
-  Fix barrier API documentation
-  On Apple, link Eigen against Accelerate as a BLAS/LAPACK backend

   -  Requires ``brew install lapack`` to install LAPACKE headers

v1.3.0 (Aug 03, 2024)
---------------------

-  Separated the collision set and potential computations. This allows us to more easily add new potentials in the future. This will require updating calls to ``compute_potential_*``. See the `tutorial <https://ipctk.xyz/tutorial/getting_started.html>`__ for details.
-  Add a ``Barrier`` class to enable dynamic selection of barrier function.
-  Add a ``NarrowPhaseCCD`` class to enable dynamic selection of narrow-phase CCD method.

.. _details-4:

Details
~~~~~~~

-  Refactor potentials in `#83 <https://github.com/ipc-sim/ipc-toolkit/pull/83>`__

   -  Replace "constraint" and "contact" names with "collision"
   -  Removed ``compute_potential_*`` from ``Collision`` and ``Collisions``
   -  Add a new class hierarchy ``Potential`` which represents computing the sum of individual potentials per collision

      -  Implement the barrier potential and friction dissipative potential as Potentials: ``BarrierPotential`` and ``FrictionPotential``

   -  Now, ``Collisions`` serve solely as the set of active collisions
   -  Add the distance mollifier to all collisions with a ``is_mollified()`` function

      -  The default mollifier is $m(x) = 1$, and only ``EdgeEdgeCollision`` overrides this

   -  Remove versions of ``compute_distance`` and ``ccd`` from ``CollisionStencil`` which take the full mesh as input

      -  Instead, expose the versions that take the collision stencil's vertex positions directly

-  Polymorphic barrier in `#84 <https://github.com/ipc-sim/ipc-toolkit/pull/84>`__

   -  Make the barrier function an object so it can be changed at runtime
   -  Add a virtual class ``Barrier`` as an interface for generic barriers
   -  Add ``ClampedLogBarrier`` class which implements the smoothly clamped log barrier functions from [Li et al. 2020]
   -  ``barrier_gradient`` and ``barrier_hessian`` renamed to ``barrier_first_derivative`` and ``barrier_second_derivative`` respectively
   -  Co-authored by `@arvigj <https://github.com/arvigj>`__

-  Update compiler warnings in `#88 <https://github.com/ipc-sim/ipc-toolkit/pull/88>`__

   -  Add ``-Werror=enum-conversion`` and ``-Wfloat-conversion``

-  Fix 🐛 hash when not using Abseil in `#90 <https://github.com/ipc-sim/ipc-toolkit/pull/90>`__
-  Clean-up SpatialHash Broad Phase in `#91 <https://github.com/ipc-sim/ipc-toolkit/pull/91>`__

   -  Replace ``IPC_TOOLKIT_WITH_CORRECT_CCD`` with ``IPC_TOOLKIT_WITH_INEXACT_CCD``

      -  Always include Tight Inclusion CCD because it is used by Nonlinear CCD

   -  Add support for face-face collision (not used anywhere but for completeness and future-proof nonlinear face-face)
   -  Add and use generic functions to ``SpatialHash``
   -  Replace ``camelCase`` with ``snake_case`` in ``SpatialHash`` and ``HashGrid``

-  Update Scalable CCD in `#92 <https://github.com/ipc-sim/ipc-toolkit/pull/92>`__

   -  Updated Scalable CCD (i.e., Sweep and Tiniest Queue and CUDA Tight Inclusion CCD) to the `unified repository <https://github.com/Continuous-Collision-Detection/Scalable-CCD>`__ with support for generic collision pairs.
   -  Renamed ``SWEEP_AND_TINIEST_QUEUE`` to ``SWEEP_AND_PRUNE`` to reflect that it is a standard implementation of the sweep and prune algorithm (see, e.g., "Real-Time Collision Detection" [Ericson 2004])
   -  Renamed ``SWEEP_AND_TINIEST_QUEUE_GPU`` to ``SWEEP_AND_TINIEST_QUEUE`` to reflect that it is the only existing implementation of the sweep and tiniest queue algorithm

-  Mark single argument constructors explicit in `#93 <https://github.com/ipc-sim/ipc-toolkit/pull/93>`__

   -  Mark ``BarrierPotential(double)`` and ``FrictionPotential(double)`` as explicit constructors to avoid implicit conversions from double

-  Fix compatibility with the latest Eigen by `@teseoch <https://github.com/teseoch>`__ in `#94 <https://github.com/ipc-sim/ipc-toolkit/pull/94>`__
-  Add clang-format check action in `#96 <https://github.com/ipc-sim/ipc-toolkit/pull/96>`__
-  Add new project PSD option by `@Huangzizhou <https://github.com/Huangzizhou>`__ in `#95 <https://github.com/ipc-sim/ipc-toolkit/pull/95>`__

   -  Instead of clamping the negative eigenvalues to zero, add the option to flips the sign of negative eigenvalues according to `[Chen et al. 2024] <https://github.com/honglin-c/abs-psd>`__

-  Fix 🐛 Python bindings for Numpy 2 in `#100 <https://github.com/ipc-sim/ipc-toolkit/pull/100>`__
-  Make faces an optional parameter to ``CollisionMesh`` in `#105 <https://github.com/ipc-sim/ipc-toolkit/pull/105>`__
-  Fix Python documentation by `@rc <https://github.com/rc>`_  in `#108 <https://github.com/ipc-sim/ipc-toolkit/pull/108>`__
-  Polymorphic CCD in `#110 <https://github.com/ipc-sim/ipc-toolkit/pull/110>`__

   -  Add narrow phase CCD parent class; pass CCD object to choose method
   -  Replace ``tolerance`` and ``max_iterations`` parameters with ``const NarrowPhaseCCD& narrow_phase_ccd`` parameter
   -  ``NarrowPhaseCCD`` is a virtual class containing the CCD methods for point-point, point-edge, edge-edge, and point-triangle CCD
   -  ``NarrowPhaseCCD`` is implemented by ``InexactCCD``, ``TightInclusionCCD``, and ``AdditiveCCD`` classes
   -  **[Breaking]** The optional parameter order to ``is_step_collision_free`` and ``compute_collision_free_stepsize`` is changed from

      .. code:: cpp

         const BroadPhaseMethod broad_phase_method = DEFAULT_BROAD_PHASE_METHOD,
         const double min_distance = 0.0,
         const double tolerance = DEFAULT_CCD_TOLERANCE,
         const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS);

      to

      .. code:: cpp

         const double min_distance = 0.0,
         const BroadPhaseMethod broad_phase_method = DEFAULT_BROAD_PHASE_METHOD,
         const NarrowPhaseCCD& narrow_phase_ccd = DEFAULT_NARROW_PHASE_CCD);

   -  The inexact floating-point CCD can be enabled beside the Tight Inclusion CCD rather than replacing it

v1.2.1 (Jul 12, 2024)
---------------------

Bug fixes |:bug:| :

- Update Pybind11 to support Numpy 2.0. Fixes segmentation fault as described in `#102 <https://github.com/ipc-sim/ipc-toolkit/issues/102>`__.

v1.2.0 (Dec 11, 2023)
---------------------

Various new features |:rocket:| and some bug fixes |:bug:|.

-  Implement the improved max approximator as described in `[Li et al. 2023] <https://arxiv.org/abs/2307.15908>`__
-  Add a port of the Additive CCD method from `[Li et al. 2021] <https://ipc-sim.github.io/C-IPC/>`__
-  Add a generic implementation of the nonlinear CCD (of linear geometry) algorithm from `[Ferguson et al. 2021] <https://ipc-sim.github.io/rigid-ipc/>`__
-  Add missing codimensional collision support (point-point and point-edge)

.. _details-3:

Details
~~~~~~~

* Update website URL to `ipctk.xyz <https://ipctk.xyz>`__ in `#54 <https://github.com/ipc-sim/ipc-toolkit/pull/54>`__
* Simplify tangential basis Jacobian calculation thanks to `@halehOssadat <https://github.com/halehOssadat>`__ and `@jpanetta <https://github.com/jpanetta>`__ in `#56 <https://github.com/ipc-sim/ipc-toolkit/pull/56>`__
* Update ``FindSIMD.cmake`` to now add support for Neon (Arm/Apple Silicon SIMD instruction set) in `#58 <https://github.com/ipc-sim/ipc-toolkit/pull/58>`__

  * Credit: ``FindSIMD.cmake`` from `Project CHRONO <https://github.com/projectchrono/chrono>`__ under `BSD 3-Clause “New” or “Revised” License <https://github.com/projectchrono/chrono/blob/main/LICENSE>`__.

* Improve the max approximator used (i.e., sum over constraints) as described in `[Li et al. 2023] <https://arxiv.org/abs/2307.15908>`__ in `#55 <https://github.com/ipc-sim/ipc-toolkit/pull/55>`__

  * Add a ``dtype`` to EE collisions to keep track of the distance type for mollified constraints
  * Initialize mesh adjacencies by default
  * Use edge length as the area weighting for codimensional edges

* Improve documentation and tutorials in `#61 <https://github.com/ipc-sim/ipc-toolkit/pull/61>`__

  * Add documentation describing the convergent formulation
  * Add documentation describing the constraint offset/minimum distance
  * Add documentation for broad- and narrow-phase CCD
  * Add documentation for High-Order IPC
  * Also, renames ``CollisionConstraint::minimum_distance`` to ``CollisionConstraint::dmin``

* Add a port of the Additive CCD method from `[Li et al. 2021] <https://ipc-sim.github.io/C-IPC/>`__ in `#62 <https://github.com/ipc-sim/ipc-toolkit/pull/62>`__

  * This is a modified version of the `original open-source implementation <https://github.com/ipc-sim/Codim-IPC>`__ which is under the `Appache-2.0 License <https://github.com/ipc-sim/Codim-IPC/blob/main/LICENSE>`__.
  * Modifications: remove broad phase functions, refactor code to use a single implementation of the ``additive_ccd`` algorithm, utilize our distance function rather than porting the Codim-IPC versions, return ``true`` if the initial distance is less than the minimum distance, and add an explicit ``tmax`` parameter rather than relying on the initial value of ``toi``.
  * This is mostly for reference comparison and it is not integrated into the full code. This also includes the ability to pull the sample CCD queries and run them in a unit-test (requires GMP).
  * This adds missing feature mentioned in `#63 <https://github.com/ipc-sim/ipc-toolkit/discussions/63>`__

* Add Codecov to get a report of unit test code coverage in `#64 <https://github.com/ipc-sim/ipc-toolkit/pull/64>`__

  * Add more tests to improve code coverage and fix small bugs in `#65 <https://github.com/ipc-sim/ipc-toolkit/pull/65>`__

* Fix the symmetric matrix assertion in ``project_to_psd`` and ``project_to_pd`` in `#67 <https://github.com/ipc-sim/ipc-toolkit/pull/67>`__
* Handle codim. point-point collisions in `#66 <https://github.com/ipc-sim/ipc-toolkit/pull/66>`__

  * This adds missing feature as discussed in `#63 <https://github.com/ipc-sim/ipc-toolkit/discussions/63>`__

* Add tests of Python bindings using `nose2 <https://docs.nose2.io/en/latest/>`__ in `#69 <https://github.com/ipc-sim/ipc-toolkit/pull/69>`__
* In CCD, check the initial distance when no motion occurs in `#71 <https://github.com/ipc-sim/ipc-toolkit/pull/71>`__
* Add a generic implementation of the nonlinear CCD (of linear geometry) algorithm from `[Ferguson et al. 2021] <https://ipc-sim.github.io/rigid-ipc/>`__ in `#72 <https://github.com/ipc-sim/ipc-toolkit/pull/72>`__

  * Generic nonlinear trajectories are specified through a ``NonlinearTrajectory`` virtual class. By default the maximum distance between the trajectory and a linearized version is computed using interval arithmetic. That is

    .. math::

      \max_{t \in [0, 1]} \Vert p(\mathrm{lerp}(t_0, t_1, t)) - \mathrm{lerp}(p(t_0), p(t_1), t) \Vert_2 \\
      \leq \sup(\Vert p([t_0, t_1]) - \mathrm{lerp}(p(t_0), p(t_1), [0, 1]) \Vert_2)

    where :math:`p` is the point's position over time, :math:`\mathrm{lerp}(a, b, t) := (b - a) t + a` and :math:`\sup([a,b]):=b`. Because this can be an overly conservative approximation, users can override the ``NonlinearTrajectory::max_distance_from_linear`` function to compute the max directly in closed form, if known.
  * We perform interval arithmetic using `filib <https://github.com/zfergus/filib>`__ which has been shown to be “the only library that is correct, consistent, portable, and efficient” `[Tang et al. 2022] <https://cims.nyu.edu/gcl/papers/2022-Intervals.pdf>`__.
  * Add a nonlinear CCD tutorial to the docs in `#78 <https://github.com/ipc-sim/ipc-toolkit/pull/78>`__

* Add additional compiler warnings and resolve them to be warning-free in `#73 <https://github.com/ipc-sim/ipc-toolkit/pull/73>`__
* Add Python bindings for ``igl::predicate::segment_segment_intersect`` in `#74 <https://github.com/ipc-sim/ipc-toolkit/pull/74>`__
* Integrate `SimpleBVH <https://github.com/ipc-sim/SimpleBVH>`__ as a broad-phase method in `#75 <https://github.com/ipc-sim/ipc-toolkit/pull/75>`__
* Fix the shape derivative of mollified edge-edge contact in `#76 <https://github.com/ipc-sim/ipc-toolkit/pull/76>`__

  * Additionally, this makes the shape derivative computation object-oriented.

* Update Python bindings with recent changes and unified comments in `#77 <https://github.com/ipc-sim/ipc-toolkit/pull/77>`__
* Add support for collision between codimensional edges and points in 3D in `#79 <https://github.com/ipc-sim/ipc-toolkit/pull/79>`__

  * Implements missing features discussed in `#63 <https://github.com/ipc-sim/ipc-toolkit/discussions/63>`__.

v1.1.1 (Aug 18, 2023)
---------------------

* Logo by `@zfergus <https://github.com/zfergus>`__ in `#52 <https://github.com/ipc-sim/ipc-toolkit/pull/52>`__
* Fix vertex-vertex :cpp:`==` and :cpp:`<` functions to be order independent

  * This allows vertex-vertex constraints merge correctly

* Update Tight Inclusion CCD

v1.1.0 (Jul 25, 2023)
---------------------

Large refactoring to make the code more object-oriented rather than passing objects to functions. Other changes include the friction potential now being a function of velocity, bug fixes, and a new tutorial.

.. _details-2:

Details
~~~~~~~

* Large Refactor in `#25 <https://github.com/ipc-sim/ipc-toolkit/pull/25>`__

  * :cpp:`construct_collision_candidates(..., candidates)` → :cpp:`candidates.build(...)`
  * :cpp:`is_step_collision_free(candidates, ...)` → :cpp:`candidates.is_step_collision_free(...)`
  * :cpp:`compute_collision_free_stepsize(candidates, ...)` → :cpp:`candidates.compute_collision_free_stepsize(...)`
  * :cpp:`compute_barrier_potential*(constraints, ...)` → :cpp:`constraints.compute_potential*(...)`
  * :cpp:`compute_shape_derivative(constraints, ...)` → :cpp:`constraints.compute_shape_derivative(...)`
  * :cpp:`compute_minimum_distance(constraints, ...)` → :cpp:`constraints.compute_minimum_distance(...)`
  * :cpp:`construct_friction_constraint_set(..., friction_constraints)` → :cpp:`friction_constraints.build(...)`
  * :cpp:`compute_friction_*(..., friction_constraints, ...)` → :cpp:`friction_constraints.compute_*(...)`
  * Generic :cpp:`CollisionStencil` parent class to :cpp:`Candidates`, :cpp:`CollisionConstraints`, and :cpp:`FrictionConstraints`.
  * Renamed :cpp:`Constraints` to :cpp:`CollisionConstraints`
  * Replaced single letter variable names :cpp:`V`, :cpp:`E`, :cpp:`F` with :cpp:`vertices`/:cpp:`positions`, :cpp:`edges`, :cpp:`faces`
  * Renamed ``*_index`` → ``*_id``
  * Replaced :cpp:`inflation_radius = min_distance / 1.99` with :cpp:`inflation_radius = min_distance / 2` and use rounding mode to conservativly inflate AABBs
  * :cpp:`CollisionConstraints::use_convergent_formulation` and :cpp:`are_shape_derivatives_enabled` must now be accessed through getter and setter functions
  * Friction potentials are now functions of velocity. Previously :cpp:`V0` and :cpp:`V1` were passed and :cpp:`U = V1-V0`. This limited the integration scheme to implicit Euler. Upstream this means you need to multiply the potential by :math:`1/(dv/dx)` to get the correct friction force.

    * Change input :math:`\epsilon_vh` to :math:`\epsilon_v` in `#37 <https://github.com/ipc-sim/ipc-toolkit/pull/37>`__ to reflect the fact that friction is defined in terms of velocity instead of displacement now.

* Changed default :cpp:`project_hessian_to_psd` to :cpp:`false` in `#30 <https://github.com/ipc-sim/ipc-toolkit/pull/30>`__
* Update website with a tutorial (`#31 <https://github.com/ipc-sim/ipc-toolkit/pull/31>`__) and version dropdown list (`#34 <https://github.com/ipc-sim/ipc-toolkit/pull/34>`__)
* Switch from templates to using :cpp:`Eigen::Ref` in `#28 <https://github.com/ipc-sim/ipc-toolkit/pull/28>`__
* Speed up the CCD by limiting the maximum minimum distance to :cpp:`1e-4` in `#43 <https://github.com/ipc-sim/ipc-toolkit/pull/43>`__
* Fix the bug pointed out in `#41 <https://github.com/ipc-sim/ipc-toolkit/pull/41>`__ in `#42 <https://github.com/ipc-sim/ipc-toolkit/pull/42>`__. Namely, to get units of distance in the barrier we should divide the original function by :math:`\hat{d}\cdot(\hat{d} + 2d_{\min})^2` when using distance squared. Before it was being divided by :math:`2d_{\min} \hat{d} + \hat{d}^2`.
* Fix build for IPC_TOOLKIT_WITH_CORRECT_CCD=OFF in `#44 <https://github.com/ipc-sim/ipc-toolkit/pull/44>`__
* Switched from FetchContent to CPM in `#48 <https://github.com/ipc-sim/ipc-toolkit/pull/48>`__. This provides better caching between builds. Additionally, made robin-map and Abseil optional dependencies.
* Add the CFL-Inspired Culling of CCD as described in Section 3 of the Technical Supplement to IPC in `#50 <https://github.com/ipc-sim/ipc-toolkit/pull/50>`__

v1.0.0 (Feb 21, 2023)
---------------------

This is the first official release. |:rocket:|

This is a stable release of the toolkit prior to refactoring the code and making updates to the API.

.. _details-1:

Details
~~~~~~~

* Added a minimum distance optional parameter to all CCD functions (:cpp:`const double min_distance = 0.0`) in `#22 <https://github.com/ipc-sim/ipc-toolkit/pull/22>`__. This is placed as the first optional argument which can break calling code if optional parameters were previously used.
* Added :cpp:`CollisionMesh` in `#7 <https://github.com/ipc-sim/ipc-toolkit/pull/7>`__ to wrap up face and edges into a single data structure.

  * Removes Support for ignoring internal vertices. Instead, users should use the CollisionMesh to map from the full mesh to the surface mesh.
  * This also includes a :cpp:`to_full_dof` function that can map the reduced gradient/hessian to the full mesh's DOF.

Pre-v1.0.0
----------

2021-10-05 (`9e2cc2a <https://github.com/ipc-sim/ipc-toolkit/commit/574f7577daa5e0b51bf5baf20998994b8371216e>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Added
^^^^^

* Added implicits source folder to organize point-plane collisions

.. _e2cc2a-1:

2021-09-05 (`9e2cc2a <https://github.com/ipc-sim/ipc-toolkit/commit/9e22cc2a5f7e7ca048a579f2c94d2241782ecf17>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-1:

Added
^^^^^

* Added support for point vs. (static) analytical plane contact

2021-08-21 (`acf2a80 <https://github.com/ipc-sim/ipc-toolkit/commit/acf2a80544ebe27dc5e440602a3a89243e575e8a>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Changed
^^^^^^^

* Changed CMake target name to ``ipc::toolkit``

2021-07-26 (`1479aae <https://github.com/ipc-sim/ipc-toolkit/commit/1479aaea958daaa4e963529493e4169dc7757913>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _changed-1:

Changed
^^^^^^^

* Updated the CMake system to use modern ``FetchContent`` to download externals

2021-07-22 (`e24c76d <https://github.com/ipc-sim/ipc-toolkit/commit/e24c76ddc818fb9efc4d522ef72a581a15abf751>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fixed
^^^^^

* Updated CCD strategy when using Tight Inclusion to only perform :cpp:`no_zero_toi=true` when there is no minimum distance

2021-07-17 (`a20f7a2 <https://github.com/ipc-sim/ipc-toolkit/commit/a20f7a2dfea5a04c67ef71d0cd523f69391f2f54>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-2:

Added
^^^^^

* Added :cpp:`detect_edge_face_collision_candidates_brute_force` for 3D intersection broad-phase
* Added ability to save an obj of collision candidates
* Added tests for has_intersection (all pass after fixes)

.. _fixed-1:

Fixed
^^^^^

* Fixed possible numerical rounding problems in HashGrid :cpp:`AABB::are_overlapping`
* Fixed HashGrid's function for getting edge-face intersection candidates

2021-07-15 (`7301b42 <https://github.com/ipc-sim/ipc-toolkit/commit/7301b422a9b9a90c76d9e7abf2f9127bf6d0dbd6>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _fixed-2:

Fixed
^^^^^

* Use :cpp:`ignore_codimensional_vertices` in the brute force broad-phase method
* Fixed AABB inflation in brute force and SpatialHash methods

2021-07-08 (`86ae4e5 <https://github.com/ipc-sim/ipc-toolkit/commit/86ae4e5f87eb2c65585920ad3ca0bbb3b57702f6>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _changed-2:

Changed
^^^^^^^

* Replaced vertex group ids with more powerful can_collide function. By default everything can collide with everything (same as before)
* Reordered parameters in :cpp:`construct_constraint_set()`, :cpp:`is_collision_free()`, and :cpp:`compute_collision_free_stepsize()`
* :cpp:`update_barrier_stiffness` now requires the :cpp:`constraint_set` rather than building it
* :cpp:`update_barrier_stiffness` dropped dhat parameter

.. _fixed-3:

Fixed
^^^^^

* SpatialHash for 2D

Removed
^^^^^^^

* Verison of :cpp:`initial_barrier_stiffness` that computes the constraint set and barrier gradient because there are a lot of parameters to these functions

2021-07-05 (`4d16954 <https://github.com/ipc-sim/ipc-toolkit/commit/4d16954012570b3a15346b99b5aedea77266fe86>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _changed-3:

Changed
^^^^^^^

* Renamed directory ``src/spatial_hash/`` → ``src/broad_phase/``
* Renamed files ``src/ccd/broad_phase.*`` → ``src/ccd/aabb.*``

2021-07-05 (`b3808e1 <https://github.com/ipc-sim/ipc-toolkit/commit/b3808e15bdbaba9a6efd4b731db3070e85bcc4b7>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-3:

Added
^^^^^

* Select the broad-phase method for CCD and distance constraints

  * Methods: :cpp:`HASH_GRID`, :cpp:`SPATIAL_HASH`, :cpp:`BRUTE_FORCE`

* CCD parameters for Tight Inclusion's tolerance and maximum iterations

.. _changed-4:

Changed
^^^^^^^

* :cpp:`ignore_codimensional_vertices` to :cpp:`false` by default
* CMake option ``TIGHT_INCLUSION_WITH_NO_ZERO_TOI=ON`` as default

2021-06-18 (`aa59aeb <https://github.com/ipc-sim/ipc-toolkit/commit/aa59aeb0634af981a8f1cfbb6d2ff2b76a04d610>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _changed-5:

Changed
^^^^^^^

* :cpp:`construct_friction_constraint_set` now clears the given :cpp:`friction_constraint_set`

2021-05-18 (`245b13b <https://github.com/ipc-sim/ipc-toolkit/commit/245b13bcc5e99ed52850ae865aaa0ad4e71a43a8>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _changed-6:

Changed
^^^^^^^

* Use TightInclusion degenerate edge-edge for point-point and point-edge CCD

2021-05-11 (`5c34dcd <https://github.com/ipc-sim/ipc-toolkit/commit/5c34dcdf226d46ada962204585fa386eb9b67859>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _changed-7:

Changed
^^^^^^^

* :cpp:`char*` exceptions to :cpp:`std::exceptions`

2021-05-06 (`24056cc <https://github.com/ipc-sim/ipc-toolkit/commit/24056ccb2ca0a03bdef8141bc5011c41547f06b5>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _changed-8:

Changed
^^^^^^^

* Gave :cpp:`dhat_epsilon_scale` a default value of :cpp:`1e-9` in :cpp:`update_barrier_stiffness`
* :warning: Changed order of parameters to :cpp:`update_barrier_stiffness`

  * Flipped :cpp:`bbox_diagonal` and :cpp:`dhat_epsilon_scale`

2021-05-06 (`81d65f3 <https://github.com/ipc-sim/ipc-toolkit/commit/81d65f32e479fea32d0acc29c8a7a532fa55518b>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _fixed-4:

Fixed
^^^^^

* Bug in output min distance of :cpp:`update_barrier_stiffness`

2021-05-04 (`59ec167 <https://github.com/ipc-sim/ipc-toolkit/commit/59ec167b85eaf56095a2d0333bdd96146d658ebf>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _changed-9:

Changed
^^^^^^^

* Moved eigen_ext functions into ipc namespace
* Renamed max size matrices with ``Max``

  * ``Eigen::VectorX([0-9])`` → ``ipc::VectorMax$1``
  * ``Eigen::MatrixXX([0-9])`` → ``ipc::VectorMax$1``
  * ``Eigen::ArrayMax([0-9])`` → ``ipc::ArrayMax$1``

2021-05-03 (`664d65f <https://github.com/ipc-sim/ipc-toolkit/commit/664d65fd70dbd350b6bfe5f8a311a89ff4fef3bd>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-4:

Added
^^^^^

* Added utility function to check for edge-edge intersection in 2D and edge-triangle intersection in 3D.
* Optionally: use GMP for exact edge-triangle intersection checks

2021-05-03 (`9b4ebfc <https://github.com/ipc-sim/ipc-toolkit/commit/9b4ebfc0f458645cf33eeebf8211607f45ad9cb4>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-5:

Added
^^^^^

* voxel_size_heuristic.cpp which suggests a good voxel size for the :cpp:`SpatialHash` and :cpp:`HashGrid`

.. _changed-10:

Changed
^^^^^^^

* Changed HashGrid voxel size to be the average edge length not considering displacement length. This results in better performance, but can result in large memory usage.

2021-04-29 (`293d0ad <https://github.com/ipc-sim/ipc-toolkit/commit/293d0ad992c01df561e25c286043c9ae9b901ff0>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-6:

Added
^^^^^

* Added TBB parallel loops to the main function (:cpp:`compute_potential`, :cpp:`compute_friction_potential`, :cpp:`compute_collision_free_stepsize`, etc.)
* Added function :cpp:`addVerticesFromEdges` that adds the vertices connected to edges in parallel and avoids duplicates

.. _changed-11:

Changed
^^^^^^^

* Changed the HashGrid to use :cpp:`ArrayMax3` over :cpp:`VectorX3` to simplify the code

.. _fixed-5:

Fixed
^^^^^

* Fixed some parameters that were not by reference

2021-04-21 (`c8a6d5 <https://github.com/ipc-sim/ipc-toolkit/commit/c8a6d56823793e7be5e89238c3793e25bc45ffa0>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-7:

Added
^^^^^

* Added the SpatialHash from the original IPC code base with some modification to get all candidates in parallel

  * Benchmark results indicate this SpatialHash is faster than the HashGrid with multithreading
  * TODO: Improve HashGrid or fully integrate SpatialHash into ipc.hpp

2021-02-11 (`9c7493 <https://github.com/ipc-sim/ipc-toolkit/commit/9c74938fefa691db6b79c73489c8c661638019c6>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _changed-12:

Changed
^^^^^^^

* Switched to the correct (conservative) CCD of :cite:t:`Wang2021TightInclusion`

  * Can select Etienne Vouga's CCD in the CMake (see README.md)

2021-02-01 (`b510253 <https://github.com/ipc-sim/ipc-toolkit/commit/b51025310223b487e7c39858265d8d5c3e8b1e8a>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-8:

Added
^^^^^

* Added minimum seperation distance (thickness) to distance constraints

  * Based on "Codimensional Incremental Potential Contact" :cite:p:`Li2021CIPC`

2021-02-01 (`a395175 <https://github.com/ipc-sim/ipc-toolkit/commit/a3951750ca5f167ab1d546ae1dadd87d0a9e2497>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-9:

Added
^^^^^

* Added 2D friction model based on the 3D formulation.

  * TODO: Test this further

2021-01-12 (`deee6d0 <https://github.com/ipc-sim/ipc-toolkit/commit/deee6d0f9802910c5565f800492f9a995e65cf7e>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-10:

Added
^^^^^

* Added and optional parameter :cpp:`F2E` to :cpp:`construct_constraint_set()`. This is similar to :cpp:`F` (which maps faces to vertices), but maps faces to edges. This is optional, but recommended for better performance. If not provided a simple linear search will be done per face edge!

  * TODO: Add a function to compute this mapping.

.. _deee6d0-1:

2021-01-09 (`deee6d0 <https://github.com/ipc-sim/ipc-toolkit/commit/deee6d0f9802910c5565f800492f9a995e65cf7e>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _changed-13:

Changed
^^^^^^^

* Replaced VectorXd and MatrixXd with static size versions for local gradient and hessians

2020-11-20 (`93143ad <https://github.com/ipc-sim/ipc-toolkit/commit/93143ad9b31030cde7324a83354268021e1cb9da>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _changed-14:

Changed
^^^^^^^

* Removed TBB parallelization form the hash grid because we get better performance without it.

  * TODO: Improve parallelization in the hash grid or switch to the original IPC spatial hash

2020-11-06 (`4553509 <https://github.com/ipc-sim/ipc-toolkit/commit/4553509fe6a4e6b78c041018cd6db3fdf23b4730>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _fixed-6:

Fixed
^^^^^

* Fixed multiplicity for point-triangle distance computation to avoid duplicate point-point and point-edge pairs.

2020-10-22 (`51f4903 <https://github.com/ipc-sim/ipc-toolkit/commit/51f49030dbeec15a6a7544826f5531811a779402>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _fixed-7:

Fixed
^^^^^

* Projection of the hessian to PSD. This was completely broken as the projected matrix was never used.

2020-10-22 (`9be6c0f <https://github.com/ipc-sim/ipc-toolkit/commit/9be6c0f7e2534e426e3f09f4c547406d50d5cf9c>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _fixed-8:

Fixed
^^^^^

* Mollification of EE constraints that have a distance type of PP or PE
* If there is no mollification needed then the PP and PE constraints are stored with multiplicity
* Set the parallel EE friction constraint threshold to eps_x like in IPC

  * This avoid needing the mollification for the normal force and these forces are small anyways

2020-10-10 (`cb8b53f <https://github.com/ipc-sim/ipc-toolkit/commit/cb8b53fb098598ba5e8c95d4bdb4730e8df9382e>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _fixed-9:

Fixed
^^^^^

* Assertions in :cpp:`compute_collision_free_stepsize`

2020-10-10 (`4a5f84f <https://github.com/ipc-sim/ipc-toolkit/commit/4a5f84f1177bdae1a265dc15a84603bbc389936d>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _fixed-10:

Fixed
^^^^^

* Point-triangle distance type by replacing it with the one used in the original IPC code

2020-10-10 (`1d51a61 <https://github.com/ipc-sim/ipc-toolkit/commit/1d51a61d60bb25e08c9937285ff9e44459a2223f>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-11:

Added
^^^^^

* Boolean parameter in :cpp:`compute_friction_potential_hessian` that controls if the hessian is projected to PSD

2020-10-09 (`b737fb0 <https://github.com/ipc-sim/ipc-toolkit/commit/b737fb0e708eac5a7775766f162a5d2067db2fa4>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-12:

Added
^^^^^

* Parameter for vertex group IDs to exclude some collisions (e.g., self collisions)

2020-10-08 (`6ee60ae <https://github.com/ipc-sim/ipc-toolkit/commit/6ee60aeaef6d7f88013ee2ee3d544e7403282527>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-13:

Added
^^^^^

* Second version of :cpp:`update_barrier_stiffness()` that takes an already computed minimum distance and world bounding box diagonal

2020-10-08 (`cc3947d <https://github.com/ipc-sim/ipc-toolkit/commit/cc3947d48bc069488f6a773424e30fe67eb4b5f1>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-14:

Added
^^^^^

* Second version of :cpp:`initial_barrier_stiffness()` that takes an already computed barrier gradient
* Assertions on :cpp:`initial_barrier_stiffness()` input

  * :cpp:`average_mass > 0 && min_barrier_stiffness_scale > 0`

.. _changed-15:

Changed
^^^^^^^

* Fixed typo in :cpp:`initial_barrier_stiffness()` name (was :cpp:`intial_barrier_stiffness()`)

.. _section-1:

2020-10-07 (`5582582 <https://github.com/ipc-sim/ipc-toolkit/commit/5582582bc2f54464bfcee4ba0ec2b7e6975f596f>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-15:

Added
^^^^^

* :cpp:`FrictionConstraint` structures to store friction information (i.e., tangent basis, normal force magnitude, closest points, and coefficient of friction)
* Unit test that compares the original IPC code's friction components with the toolkit's

.. _changed-16:

Changed
^^^^^^^

* :cpp:`compute_friction_bases()` is now :cpp:`construct_friction_constraint_set()`

  * It now takes the coefficient of friction (:cpp:`mu`)
  * It now puts all information inside of the :cpp:`FrictionConstraints` (:cpp:`friction_constraint_set`)

2020-10-06 (`b48ba0e <https://github.com/ipc-sim/ipc-toolkit/commit/b48ba0ec9d60754e7670e28fd1987b0c78cd809f>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _changed-17:

Changed
^^^^^^^

* During :cpp:`construct_constraint_set()` the constraints are added based on distance type

  * Duplicate vertex-vertex and edge-vertex constraints are handled by a multiplicity multiplier
  * Edge-edge constraints are always line-line distances
  * Point-triangle constraints are always point-plane distances

2020-10-05 (`9a4576b <https://github.com/ipc-sim/ipc-toolkit/commit/9a4576b209302c79296593ac213ed8ce85510f3b>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _fixed-11:

Fixed
^^^^^

* Fixed a bug in the point-triangle closest points and tangent basis computed in :cpp:`compute_friction_bases()`
* Fixed a bug in :cpp:`edge_edge_tangent_basis()` used to compute the tangent basis for friction

2020-09-19 (`31a37e0 <https://github.com/ipc-sim/ipc-toolkit/commit/31a37e04abc9ecec325e00be97fd42b89c895b45>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-16:

Added
^^^^^

* spdlog for logging information

2020-09-19 (`acb7664 <https://github.com/ipc-sim/ipc-toolkit/commit/acb7664792982685f6de28468ba126f5e531834f>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _changed-18:

Changed
^^^^^^^

* Headers are now include with the prefix ``ipc/``

  * E.g., :cpp:`#include <ipc.hpp>` → :cpp:`#include <ipc/ipc.hpp>`

2020-09-04 (`7dd2ab7 <https://github.com/ipc-sim/ipc-toolkit/commit/7dd2ab7a255ffd23ccdfe5aee08bca6a142f75a7>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _added-17:

Added
^^^^^

* Collision constraint to store distance constraint pairs

  * :cpp:`EdgeEdgeConstraint` stores the edge-edge mollifier threshold (:cpp:`eps_x`)

.. _changed-19:

Changed
^^^^^^^

* Input parameter :cpp:`dhat_squared` is now :cpp:`dhat` (i.e., non-squared value)
* Input parameter :cpp:`epsv_times_h_squared` is now :cpp:`epsv_times_h` (i.e., non-squared value)
* :cpp:`Constraints` replaced :cpp:`Candidates`
* :cpp:`construct_constraint_set()` now takes the rest vertex position (:cpp:`V_rest`)
* :cpp:`compute_barrier_potential*()` no longer take the rest vertex position
