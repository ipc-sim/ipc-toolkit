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
2. Use the ``TrustRegion`` helper class to filter the step and compute the per-vertex safe distances.

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

Trust Region Logic
~~~~~~~~~~~~~~~~~~

From the above pseudocode, the ``TrustRegion`` helper class implements the logic for each piece.

Start by constructing a ``TrustRegion`` instance with the desired parameters:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            #include <ipc/ogc/trust_region.hpp>

            ipc::ogc::TrustRegion trust_region(dhat);
            // Optionally set parameters like:
            trust_region.relaxed_radius_scaling = 0.9; // 2γₚ in the paper
            trust_region.update_threshold = 0.01; // γₑ in the paper

    .. md-tab-item:: Python

        .. code-block:: python

            from ipctk import ogc

            trust_region = ogc.TrustRegion(dhat)
            # Optionally set parameters like:
            trust_region.relaxed_radius_scaling = 0.9  # 2γₚ in the paper
            trust_region.update_threshold = 0.01  # γₑ in the paper


Next step 2b in the pseudocode is to compute the conservative step for each vertex.
This is done via ``TrustRegion::compute_conservative_step``,
