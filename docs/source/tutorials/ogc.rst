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

For more details on the OGC algorithm, please refer to :cite:`Chen2025Offset`
and their associated video:

.. youtube:: xxyniqSLJik
   :width: 100%

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

            // e. Re-check for Collision Detection
            // If too many vertices moved beyond their conservative bounds,
            // a new, full collision detection is needed in the next step
            if num_vertices_exceed_bound ≥ threshold (γₖ):
                collision_detection_required = true

            // f. Check for Convergence
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
            trust_region.update_threshold = 0.01;      // γₑ in the paper

    .. md-tab-item:: Python

        .. code-block:: python

            from ipctk import ogc

            trust_region = ogc.TrustRegion(dhat)
            # Optionally set parameters like:
            trust_region.relaxed_radius_scaling = 0.9  # 2γₚ in the paper
            trust_region.update_threshold = 0.01       # γₑ in the paper


Step 1 & 2c: Initialization and Prediction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first step of the algorithm (Step 1) and the initial collision detection update (Step 2a) are handled by ``warm_start_time_step``. This function initializes the trust region around the current positions and moves the vertices towards their predicted positions (``pred_x``) while respecting the conservative bounds.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            // Initialize trust region and move x towards pred_x
            trust_region.warm_start_time_step(
                collision_mesh, x, pred_x, collisions, dhat);

    .. md-tab-item:: Python

        .. code-block:: python

            # Initialize trust region and move x towards pred_x
            trust_region.warm_start_time_step(
                collision_mesh, x, pred_x, collisions, dhat)

This function effectively performs the initial "Apply Initial Guess" and "Compute Conservative Step" logic.

Step 2d: Conservative Bound Check (Filter Step)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Inside your solver loop (e.g., Newton's method or Vertex Block Descent), you must ensure that the computed search direction or step (``dx``) does not violate the trust region bounds. This corresponds to Step 2d in the pseudocode.

Use the ``filter_step`` method to truncate any vertex displacements that exceed their trust region radius.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            // dx is the proposed step (e.g., from the linear solver)
            trust_region.filter_step(collision_mesh, x, dx);

            // Now dx is safe to apply: x += dx

    .. md-tab-item:: Python

        .. code-block:: python

            # dx is the proposed step
            trust_region.filter_step(collision_mesh, x, dx)

            # Now dx is safe to apply: x += dx

Step 2a & 2e: Re-check for Collision Detection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If ``filter_step`` detects that a significant number of vertices (determined by ``update_threshold``) were restricted by the trust region, it flags that a collision update is needed. You should check this at the beginning of each solver iteration using ``update_if_needed``.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            // At the start of the solver loop:
            trust_region.update_if_needed(collision_mesh, x, collisions, dhat);

    .. md-tab-item:: Python

        .. code-block:: python

            # At the start of the solver loop:
            trust_region.update_if_needed(collision_mesh, x, collisions, dhat)

Full Optimization Loop Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Putting it all together, a single simulation step using OGC looks like this:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            // 1. Initial setup
            ipc::NormalCollisions collisions;
            collisions.set_collision_set_type(ipc::NormalCollisions::CollisionSetType::OGC);
            ipc::ogc::TrustRegion trust_region(dhat);

            // 2. Warm start (Predict & Initialize)
            // x is current position, pred_x is x^t + dt * v^t
            trust_region.warm_start_time_step(
                mesh, x, pred_x, collisions, dhat);

            // 3. Solver Loop
            for (int i = 0; i < max_iterations; ++i) {
                // Update trust region if too many vertices hit the bound in previous step
                trust_region.update_if_needed(mesh, x, collisions, dhat);

                // Compute search direction (Solver specific)
                Eigen::MatrixXd dx = compute_search_direction(x, ...);

                // Filter the step to respect OGC bounds
                trust_region.filter_step(mesh, x, dx);

                // Update positions
                x += dx;

                // Check convergence...
            }

    .. md-tab-item:: Python

        .. code-block:: python

            # 1. Initial setup
            collisions = ipctk.NormalCollisions()
            collisions.collision_set_type = ipctk.NormalCollisions.CollisionSetType.OGC
            trust_region = ipctk.ogc.TrustRegion(dhat)

            # 2. Warm start (Predict & Initialize)
            trust_region.warm_start_time_step(
                mesh, x, pred_x, collisions, dhat)

            # 3. Solver Loop
            for i in range(max_iterations):
                # Update trust region if needed
                trust_region.update_if_needed(mesh, x, collisions, dhat)

                # Compute search direction (Solver specific)
                dx = compute_search_direction(x, ...)

                # Filter the step to respect OGC bounds
                trust_region.filter_step(mesh, x, dx)

                # Update positions
                x += dx

                # Check convergence...

Algorithm Details and Adaptations
---------------------------------

While the implementation follows the core OGC algorithm (Algorithm 3 in :cite:`Chen2025Offset`), there are several practical adaptations in ``TrustRegion`` to ensure robustness, numerical stability, and compatibility with line-search-based solvers.

Step Scaling vs. Closest-Point Projection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most significant deviation from the original paper lies in how vertices are constrained when they exceed their trust region bounds (Step 2d).

* **Original Paper Approach (Projection):** The original algorithm typically projects the proposed position to the closest point within the trust region (minimizing geometric distance). While this finds the optimal position locally within the valid region, the vector from the current position to this projected point may effectively change the update direction.
* **Toolkit Implementation (Step Scaling):** The implementation in ``filter_step`` solves for a scalar :math:`\beta \in (0, 1]` that scales the original search direction :math:`\Delta x` such that the new position lies exactly on the trust region boundary.

    That is, given the trust region center :math:`c`, current position :math:`x`, and search direction :math:`\Delta x`, it solves for :math:`\beta` such that:

    .. math::

        \| x + \beta \Delta x - c \|^2 = r^2

    This value of :math:`\beta` is then used to scale the step:

    .. code-block:: c++

        // trust_region.cpp
        // Solve || x + beta * dx - c ||^2 = r^2 for beta
        // ... (quadratic formula solution) ...
        dx.row(i).array() *= beta;

    This ensures that vertices constrained by the trust region are placed exactly on the valid boundary surface, maximizing the allowed step size without violating the constraint.

**Why this matters:**
Mixing Trust Region constraints with Line Search methods can be delicate. If you use a projection method, the resulting update vector might no longer be a *descent direction* for the energy function. This can cause subsequent line searches to fail or the solver to stagnate.

By **scaling** the step instead of projecting it, the toolkit guarantees that the modified step remains parallel to the original search direction. If the solver (e.g., Newton's method) generated a valid descent direction, this method preserves that property, ensuring compatibility with line search checks that may occur after the OGC filter is applied.


Dynamic Trust Region Inflation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Step 1 (Initialization), the implementation dynamically computes the radius used for building the collision set (``trust_region_inflation_radius``). Instead of a fixed radius, it uses the maximum displacement from the predicted motion (``pred_x - x``) combined with the offset distance ``dhat``.

.. code-block:: c++

    // trust_region.cpp
    trust_region_inflation_radius = std::max(2 * dhat, dx_norm.maxCoeff());

This follows the recommendation in Section 4.1 of the paper:

    In practice, we found an :math:`r_q` of :math:`r` plus the inertial displacement magnitude to strike a good balance between query performance and bound size.

Safety Factor and Double-Sided Contact
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When computing the per-vertex safe distances (Step 2b), the implementation applies a scaling factor and explicitly divides by 2.

.. code-block:: c++

    // trust_region.cpp
    // Use < half the safe distance to account for double sided contact
    reduced_trust_region_radii *= relaxed_radius_scaling / 2;

* **Relaxed Radius Scaling:** Corresponds to the parameter :math:`2\gamma_p` in the paper (default 0.9). It shrinks the trust region slightly to provide a numerical safety margin.
* **Division by 2:** Accounts for the worst-case scenario where two primitives move towards each other; each is allowed to move only half the available distance to guarantee no intersection.

Adaptive Warm Start
~~~~~~~~~~~~~~~~~~~

The ``warm_start_time_step`` function includes an adaptive heuristic. During the initial prediction:
1.  It moves vertices to their predicted positions ``pred_x`` if they are within the trust region.
2.  If they are outside, it effectively projects them to the boundary (using the scaling method described above).
3.  **Crucially**, if the number of vertices hitting the boundary exceeds ``update_threshold`` immediately during this warm start, it triggers a second ``update()`` call.

This prevents the solver from starting with a trust region that is already too restrictive for the bulk of the mesh, potentially saving solver iterations later.