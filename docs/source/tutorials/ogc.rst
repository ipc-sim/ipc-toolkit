Offset Geometric Contact
========================

This tutorial explains how to implement the Offset Geometric Contact
(OGC) :cite:`Chen2025Offset` algorithm using the IPC Toolkit.

OGC is a contact model designed to guarantee penetration-free simulation of
objects while incurring minimal computational overhead. It works by constructing
a volumetric shape around the object by offsetting each face along its normal
direction, which ensures orthogonal contact forces and robustly allows for a
large contact radius without introducing artifacts.

OGC integrates into the optimization process by replacing the continuous
collision detection (CCD) line search with a trust region-based approach using
vertex-specific displacement bounds.

Algorithm
---------

At a high-level, the OGC algorithm proceeds as follows:

.. code-block:: javascript

    // Input: Current positions (X), previous positions (Xt), velocities (Vt),
    //        mass matrix (M), contact radius (r), query radius (r_q).
    // Output: New positions (X) for the current step.

    function SimulationStep(X, Xt, Vt, M, r, r_q):

        // 1. Initial position and prediction
        collision_detection_required = true
        X = Xt

        // Predict new positions using an explicit integration step (e.g., a simple euler step)
        Y = X + dt * Vt

        // 2. Iterative correction loop
        for i from 0 to max_iterations:

            // a. Collision detection phase (run only when required)
            if collision_detection_required:

                // i. Update collision structures with current positions
                updateCollisions(collisions, X, r, r_q)

                // ii. Store previous positions for later comparison
                X_prev = X

                // iii. Reset the flag; no collision detection needed until further notice
                collision_detection_required = false

            // b. Compute the conservative step (b) for each vertex
            for each vertex v in set vertices:
                b(v) = computeConservativeStep(collisions, v)

            // c. Simulation and correction phase
            if i is 0:
                // Apply a special initial guess/setup for the first iteration
                X = applyInitialGuess(Xt, Vt, b)

            // Apply an update step (e.g., gradient, Newton, or vertex block descent)
            // This finds the new position X by respecting contacts and forces
            X = updatePositions(X, Xt, Y, Vt, collisions)

            // d. Conservative bound check
            num_vertices_exceed_bound = 0

            // Truncate the vertex displacements to be within the conservative bound (b_v)
            for each vertex v in set vertices:
                // If the distance moved is greater than the conservative bound
                dX = X(v) - X_prev(v)
                if norm(dX) > b(v):

                    // Truncate the new position to respect the conservative bound
                    X(v) = b(v) * dX / norm(dX) + X_prev(v)

                    // Keep track of how many vertices were truncated
                    num_vertices_exceed_bound++

            // D. Re-check for Collision Detection
            // If too many vertices moved beyond their conservative bounds,
            // a new, full collision detection is needed in the next step
            if num_vertices_exceed_bound ≥ threshold (γₖ):
                collision_detection_required = true

            // E. Check for Convergence
            // If the simulation is stable and changes are small
            if evaluateConvergence(X, Xt, Vt, collisions):
                break

        return X

Essentially, OGC replaces the CCD line search with a trust region approach
that uses per-vertex displacement bounds. If too many vertices exceed their
bounds, a new collision detection is triggered to update the contact set.

Implementation
--------------

There are two main components to implementing OGC with the IPC Toolkit:

1. Set the ``CollisionSetType`` to ``OGC`` in the ``NormalCollisions`` class to
   enable OGC filtering and construction of the collision set.
2. Use the ``compute_per_vertex_safe_distances`` helper from the ``Candidates``
   class to compute the per-vertex conservative step bounds.

Constructing Candidates
~~~~~~~~~~~~~~~~~~~~~~~

Before computing the per-vertex safe distances and building the
``NormalCollisions``, you need to construct a ``Candidates`` object that
contains the broad-phase filtered potential collisions. This involves
using a spatial acceleration structure (like a BVH or grid) to find pairs of
elements (vertices, edges, faces) that are close enough to potentially collide,
given the contact radius and query radius.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            Candidates candidates;
            candidates.build(collision_mesh, vertices, dhat, dmin);

    .. md-tab-item:: Python

        .. code-block:: python

            candidates = ipctk.Candidates()
            candidates.build(collision_mesh, vertices, dhat, dmin)

Using OGC Collision Set Type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``NormalCollisions`` class supports several collision set types via
``NormalCollisions::CollisionSetType``. To build a collision set that uses the
OGC filtering and construction logic, set the collision set type to ``OGC``
before calling ``build``.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            ipc::NormalCollisions collisions;
            collisions.set_use_area_weighting(true); // optional, for convergent form
            collisions.set_collision_set_type(ipc::NormalCollisions::CollisionSetType::OGC);
            collisions.build(collision_mesh, vertices, dhat, dmin);

    .. md-tab-item:: Python

        .. code-block:: python

            collisions = ipctk.NormalCollisions()
            collisions.use_area_weighting = True  # optional, for convergent form
            collisions.collision_set_type = ipctk.NormalCollisions.CollisionSetType.OGC
            collisions.build(collision_mesh, vertices, dhat, dmin)

.. note::
    Set ``collision_set_type`` before calling ``build``; this configures how
    ``NormalCollisions`` constructs the final collision sets (vv/ev/ee/fv)
    respecting the OGC filters and any improved-max or convergent options.

Per-Vertex Safe Distances (Large Stepping)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The IPC Toolkit exposes a helper to compute how far each vertex may move
without encountering a new collision: ``Candidates::compute_per_vertex_safe_distances``.
This computes, for each vertex, a conservative maximum step length (along any
direction) that will not introduce a new collision given the current candidate
set and an inflation radius. Using these per-vertex safe distances allows an
optimizer to take larger position updates without recomputing collision detection
at every intermediate state.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            #include <ipc/candidates/candidates.hpp>
            #include <Eigen/Core>

            // candidates: ipc::Candidates already populated (broad-phase filtered)
            // collision_mesh: ipc::CollisionMesh
            // vertices: Eigen::MatrixXd (n x 3) current positions
            double inflation_radius = dhat; // same dhat used for collisions
            double min_distance = 0.0;      // numerical tolerance

            std::vector<double> per_vertex_steps = candidates.compute_per_vertex_safe_distances(
                collision_mesh, vertices, inflation_radius, min_distance);

            // per_vertex_steps[i] is a conservative distance vertex i can move
            // before a new collision may be introduced.

    .. md-tab-item:: Python

        .. code-block:: python

            # candidates: ipctk.Candidates instance already constructed
            # collision_mesh: ipctk.CollisionMesh
            # vertices: (n,3) numpy array
            inflation_radius = dhat
            min_distance = 0.0

            per_vertex_steps = candidates.compute_per_vertex_safe_distances(
                collision_mesh, vertices, inflation_radius, min_distance)

            # per_vertex_steps is a numpy array or list of floats; use it to guide
            # optimizer step sizes per-vertex.

Parameters and tuning
---------------------

- ``inflation_radius`` / ``dhat``: Choose this to match the support of your barrier
  potential. A smaller value yields ... . A larger value yields ...
- ``min_distance`` / ``epsilon``: Some OGC helpers allow a minimum distance tolerance to
  account for numerical error; choose a small positive value consistent with
  your solver tolerances.
- distance typing: When available, prefer using the same ``DistanceType`` as the
  rest of your pipeline (squared vs. euclidean) to avoid extra conversions.

Performance notes
-----------------
