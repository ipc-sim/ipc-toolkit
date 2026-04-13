.. _geometric-contact-potential-formulation:

Geometric Contact Potential
===========================

In addition to the original implementation of :cite:t:`Li2020IPC`, we also implement the Geometric Contact Potential (GCP) from :cite:t:`Huang2025GCP`.

Motivation
----------

A representative frictionless contact problem can be viewed as a geometrically constrained optimization:

.. math::
   \min_u E(u), \quad \text{subject to } g(u) \geq 0,

where :math:`u` is a deformation and :math:`g` measures the *distance to contact*. Barrier-based solvers replace this with an unconstrained problem:

.. math::
   \min_u E(u) + b(u),

where :math:`b(u)` is a barrier potential that increases to infinity as the configuration approaches contact.

All existing contact barrier potentials are defined by aggregating repulsion terms between pairs of points or elements. The key question is how to choose the strength of repulsion so that it is high when two points are close to contact, vanishes when they are far, and satisfies a number of desirable properties. In particular, one important difficulty is distinguishing *true contact* — arising from different objects or parts of the same surface moving toward each other — from points that remain close simply because they are neighbors on the undeformed rest shape.

.. figure:: /_static/img/gcp/nearby-points.svg
   :align: center
   :width: 200px
   :class: figure-on-light

   Contact potentials must distinguish points nearing contact (A and B) from nearby non-contacting points (A and C).

The IPC barrier potential :cite:p:`Li2020IPC` uses a purely distance-based criterion, which can produce *spurious repulsive forces* between nearby surface points that are not actually in contact. This is especially problematic when the barrier activation distance :math:`\hat{d}` is relatively large compared to mesh features. GCP addresses this by introducing geometric criteria — the *local minimum constraint* and the *exterior direction constraint* — that identify true contact candidates, eliminating spurious forces and decoupling :math:`\hat{d}` from the discretization.

Requirements
------------

GCP is designed to satisfy the following requirements simultaneously. No prior potential in the literature satisfies all of them:

1. **Finiteness**: The total contact potential integral is finite for any piecewise-smooth surface not in contact.
2. **Barrier**: The potential grows to infinity as the distance to contact approaches zero, guaranteeing that all configurations remain contact-free (when combined with CCD).
3. **No spurious forces**: In the undeformed (rest) configuration, both the potential and its gradient are zero, so no artificial forces arise.
4. **Localization**: The potential has a locality parameter :math:`\hat{d}` and vanishes if the distance to contact exceeds :math:`\hat{d}`.
5. **Differentiability**: The potential depends differentiably (and piecewise twice differentiably) on mesh vertex positions, enabling second-order solvers.
6. **Discretization**: The discrete version of the potential satisfies all of the above requirements exactly, not just in the limit of refinement.

Formulation
-----------

The GCP pointwise potential takes the general form:

.. math::
   \psi_\epsilon(x, y; f) := \gamma(x, y) \, p_{\epsilon(x)}(\|f(y) - f(x)\|),

where :math:`f` is the deformation map, :math:`p_{\epsilon}` is a barrier function that diverges as its argument approaches zero and vanishes for distances exceeding :math:`\epsilon(x)`, and :math:`\gamma(x, y)` is a *directional factor* that restricts the potential to a geometrically meaningful *interaction set* :math:`C(x, f)`.

The total potential is obtained by integrating (or in the discrete case, summing) over all pairs of surface points.

The formulation involves four main components:

1. **Interaction sets** :math:`C(x, f)`: the set of points that are close to being in true contact with :math:`x`.
2. **Adaptive locality** :math:`\epsilon(x)`: a per-point barrier extent.
3. **Barrier function** :math:`p_\epsilon`: the distance-based barrier.
4. **Directional factor** :math:`\gamma(x, y)`: a smooth factor supported on :math:`C(x, f)`.

Barrier Function
^^^^^^^^^^^^^^^^

The barrier function is:

.. math::
   p_\epsilon(z) := h_\epsilon(z) \, z^{-r},

where :math:`r = n - 1` (:math:`n \in \{2, 3\}` is the dimension of the scene), and :math:`h_\epsilon(z) := \tfrac{3}{2} B^3(2z/\epsilon)` is a cubic :math:`C^2` spline with support :math:`|z| \leq \epsilon`. Here :math:`B^3` is the cubic B-spline basis function:

.. math::
   B^3(v) = \begin{cases}
   \tfrac{2}{3} - v^2 + \tfrac{1}{2}|v|^3, & |v| < 1, \\
   \tfrac{1}{6}(2 - |v|)^3, & 1 \leq |v| < 2, \\
   0, & 2 \leq |v|.
   \end{cases}

Because :math:`h_\epsilon` vanishes for :math:`z \geq \epsilon` and :math:`z^{-r}` diverges as :math:`z \to 0`, the barrier satisfies both the localization and barrier requirements.

.. figure:: /_static/img/gcp/barriers.svg
   :align: center
   :width: 300px
   :class: figure-on-light

   The cubic spline :math:`h_\epsilon(z)` (dashed) and the GCP barrier :math:`p_\epsilon(z) = h_\epsilon(z)/z^{n-1}` (blue) for :math:`\epsilon = 1`. The IPC log barrier :math:`p^{\text{IPC}}` (orange) is shown for comparison.

Interaction Sets
^^^^^^^^^^^^^^^^

The key insight of GCP is defining interaction sets based on a geometric characterization of contact: *points that are close to being local minima of the distance function and for which the displacement vector points toward the surface exterior*.

For smooth surfaces, two points :math:`x` and :math:`y` on the boundary are in contact if (a) they coincide in space, :math:`f(x) = f(y)`, and (b) their normals have opposite orientation, :math:`n(x) = -n(y)`. The interaction set captures points *close* to satisfying both conditions, via two constraints:

**Local minimum constraint.** A point :math:`y` is near a local minimum of the distance from :math:`f(x)` if:

.. math::
   \Phi^m(x, y) := \|(f(y) - f(x))_+ \times n(y)\| \leq \alpha,

where :math:`(\cdot)_+` denotes normalization to a unit vector and :math:`0 < \alpha < 1`. This condition checks whether the distance direction is nearly aligned with the surface normal at :math:`y` — equivalently, whether :math:`y` is near a closest point on the surface to :math:`f(x)`.

**Exterior direction constraint.** The displacement vector :math:`f(y) - f(x)` should point toward the surface exterior at :math:`y` (i.e., :math:`f(x)` is on the outward-normal side of :math:`f(y)`):

.. math::
   \Phi^e(x, y) := -n(y) \cdot (f(y) - f(x))_+ \geq -\alpha.

This eliminates interactions between opposite sides of a thin volumetric shell.

The interaction set for smooth surfaces is then:

.. math::
   C(x, f) := \left\{ y \in \partial\Omega \;\middle|\;
   \begin{array}{l}
   \Phi^m(x, y) \leq \alpha, \quad \Phi^e(x, y) \geq -\alpha, \\
   \Phi^m(y, x) \leq \alpha, \quad \Phi^e(y, x) \geq -\alpha
   \end{array}
   \right\}.

Note that both constraints are applied symmetrically: :math:`y` must be near a local minimum as seen from :math:`x`, *and* :math:`x` must be near a local minimum as seen from :math:`y`.

.. figure:: /_static/img/gcp/demo-force.png
   :align: center
   :width: 600px
   :class: figure-on-light

   Contact forces (arrows) and potential distribution on a 2D object surface with respect to the red dot. **(A)** The IPC potential distributes spherically regardless of surface shape. **(B)** With the local minimum constraint, only surfaces near local minima/maxima of distance have high potential. **(C)** With both local minimum and exterior direction constraints, only the closer side of the volumetric object has nonzero potential.

Piecewise Smooth Surfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^

Smooth surfaces are typically approximated by piecewise-linear meshes in practice. On piecewise-smooth surfaces (including triangle meshes), normals are not uniquely defined at edges and vertices, so the constraints must be generalized.

**Local minimum constraint (generalized).** For a point :math:`y` on a piecewise-smooth surface, let :math:`I` be the set of indices of patches containing :math:`y`. For each patch :math:`\partial\Omega_i`, define tangent vectors :math:`t^k_i(y)` for :math:`k = 1, 2, 3` (two edge tangents and their angle bisector, see figure below). The local minimum condition becomes:

.. math::
   \Phi^m_{ik}(x, y) := {t^k_i(y)}^\top (f(y) - f(x))_+ \geq -\alpha_t, \quad \forall\, i \in I,\; k = 1, 2, 3,

where :math:`\alpha_t` is the tangent-direction smoothing parameter.

.. figure:: /_static/img/gcp/tangent-notation.svg
   :align: center
   :width: 600px
   :class: figure-on-light

   Tangent directions at a point :math:`f(y)` on a face :math:`\partial\Omega_i` in 3D. :math:`t_1` and :math:`t_2` are the two tangent vectors along the edge curves of the face at :math:`f(y)`, and :math:`t_3` is their angle bisector. The face is viewed from the normal direction.

**Exterior direction constraint (generalized).** At an edge point, the constraint is evaluated by projecting onto the plane perpendicular to the edge and applying a 2D inside-cone test. At a vertex, where multiple face normals :math:`n_i(y)` are available, a conservative smooth formulation is used. The per-face exterior constraint is:

.. math::
   \Phi^e_i(x, y) := -n_i(y)^\top (f(y) - f(x))_+,

and the mollified constraint ensures that :math:`\Phi^e_i(x, y) \geq 0` for *at least one* incident face :math:`i`.

.. figure:: /_static/img/gcp/normal-demo.png
   :align: center
   :width: 600px
   :class: figure-on-light

   2D (top) and 3D (bottom) normal directions for contact points. The displacement :math:`f(y) - f(x)` is shown in yellow; face normals :math:`n_i(y)` are shown in red when the exterior constraint is satisfied, and green otherwise. In simple cases such as (A) and (D), all normals satisfy the constraint, but in general only some do.

.. figure:: /_static/img/gcp/cone-filter-edge.svg
   :align: center
   :width: 600px
   :class: figure-on-light

   3D (left) and 2D (right) views of the exterior direction test for an edge point. The 2D view is the projection onto the plane perpendicular to edge :math:`e_j`. Deciding if :math:`v` points inside the cone in 3D reduces to checking if :math:`\tilde{v}` is inside the sector bounded by :math:`\tilde{e}_1` and :math:`\tilde{e}_2` in 2D.

**Contact between different element types.** On piecewise-smooth surfaces, contact can occur between any pair of element types: Face–Face, Face–Edge, Face–Vertex, Edge–Edge, Edge–Vertex, and Vertex–Vertex. Since edge and vertex contact sets have zero measure, directly integrating the potential over the surface would yield zero. GCP addresses this by treating low-dimensional elements separately with an element-measure weighting parameter :math:`L`:

.. math::
   \Psi(f) = \sum_{(g, h)} \sum_{\substack{i \in I_g,\, j \in I_h \\ G_i \cap H_j = \varnothing}} L^{4 - \dim g - \dim h} \int_{G_i} \int_{H_j} P(x, y;\, G_i, H_j) \,\mathrm{d}x\,\mathrm{d}y,

where :math:`g, h \in \{\text{Face}, \text{Edge}, \text{Vertex}\}`, and the integrals are area, line, or point evaluations depending on the element dimension.

Directional Factor :math:`\gamma`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To keep the potential differentiable while restricting it to the interaction set, GCP uses smooth mollification. The directional factor is built from two components:

**Local minimum factor** (using the smoothed Heaviside function :math:`H^\alpha`):

.. math::
   g^m(x, y) := \prod_{i \in I} \prod_{k=1}^{3} H^{\alpha_t}(\Phi^m_{ik}(x, y)).

This equals 1 when all tangent-direction constraints are satisfied, and vanishes smoothly as any constraint is violated.

**Exterior direction factor** (for vertices with multiple normals):

.. math::
   g^e(x, y) := H^1\!\left(-1 + \sum_{i \in I} H^{\beta_n}(\Phi^e_i(x, y))\right),

which equals 1 when at least one face normal satisfies the exterior constraint, and vanishes when none do. For edge points, a simpler form using the 2D cone test is used.

The combined factor is:

.. math::
   \gamma(x, y) := g^m(x, y)\, g^e(x, y)\, g^m(y, x)\, g^e(y, x).

The symmetry ensures that both points in a pair must satisfy the geometric criteria.

Smooth Heaviside Function
^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: /_static/img/gcp/heaviside.svg
   :align: center
   :width: 300px
   :class: figure-on-light

   The smooth Heaviside function :math:`H^\alpha(z)` with :math:`\alpha = \tfrac{1}{2}`.

The smooth Heaviside function :math:`H(z; \alpha, \beta) \in C^1(\mathbb{R})` used throughout the construction satisfies:

.. math::
   H(z; \alpha, \beta) = 0, \quad \forall\, z < -\alpha, \\
   H(z; \alpha, \beta) = 1, \quad \forall\, z > \beta, \\
   H'(z; \alpha, \beta) \geq 0, \quad \forall\, z.

The parameters :math:`\alpha` and :math:`\beta` control the transition region. Since the inputs to :math:`H` are dot and cross products of unit vectors, we have :math:`-1 \leq z \leq 1`, and the basic requirement is:

.. math::
   0 \leq -\beta < \alpha \leq 1.

Adaptive Barrier Localization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To satisfy the *no spurious forces* requirement, the per-point barrier extent :math:`\epsilon(x)` is set adaptively:

.. math::
   \epsilon(x) := \min\!\left(\frac{d_c(x, f_0)}{2},\; \hat{d}\right),

where :math:`d_c(x, f_0)` is the distance to the interaction set in the rest configuration and :math:`\hat{d}` is the global maximum activation distance. This ensures that the potential is exactly zero in the rest pose, eliminating spurious forces by construction.

Discretization
^^^^^^^^^^^^^^

For the barrier property to hold exactly in the discrete setting, GCP uses *closest points on element pairs as quadrature points* (rather than standard quadrature). This may reduce integration precision, but guarantees the barrier property at the discrete level. The discrete potential is:

.. math::
   \Psi^D(f) = \sum_{(g, h)} \sum_{(i, j)} L^{4 - \dim g - \dim h}\, P(x_i, y_j;\, G_i, H_j)\, A(G_i)\, A(H_j),

where :math:`(x_i, y_j)` is a pair of closest points on elements :math:`G_i` and :math:`H_j`, and :math:`A(\cdot)` denotes the element measure (area for faces, length for edges, 1 for vertices).

Distance Mollification
^^^^^^^^^^^^^^^^^^^^^^

Since closest points may not depend smoothly on vertex positions (e.g., when a closest point moves from the interior of a face to its boundary), GCP introduces mollification factors :math:`M(x, y)` that smoothly vanish near these non-smooth transitions.

.. figure:: /_static/img/gcp/smoothness-demo.svg
   :align: center
   :width: 600px
   :class: figure-on-light

   Non-smoothness of the closest point position. As the query point :math:`P` moves, the closest point :math:`Q` on an edge or face transitions from the interior to the boundary. The coordinate of :math:`Q` does not depend smoothly on :math:`P` at the boundary, motivating the mollification.

For **Edge–Vertex** contact with :math:`Q` on edge :math:`AB`:

.. math::
   M(x, y) := h_c\!\left(\frac{d(A, P)}{d(P, AB)}\right) h_c\!\left(\frac{d(B, P)}{d(P, AB)}\right),

which vanishes as the closest point approaches either endpoint of the edge.

For **Face–Vertex** contact with face :math:`ABC`:

.. math::
   M(x, y) := h_c\!\left(\frac{d(P, AB)}{d(P, ABC)}\right) h_c\!\left(\frac{d(P, AC)}{d(P, ABC)}\right) h_c\!\left(\frac{d(P, BC)}{d(P, ABC)}\right),

which vanishes as the closest point approaches the boundary of the triangle.

For **Edge–Edge** contact between edges :math:`AB` and :math:`CD`:

.. math::
   M(x, y) := h_c\!\left(\frac{d(B, CD)}{d(AB, CD)}\right) h_c\!\left(\frac{d(A, CD)}{d(AB, CD)}\right) h_c\!\left(\frac{d(C, AB)}{d(AB, CD)}\right) h_c\!\left(\frac{d(D, AB)}{d(AB, CD)}\right).

Here :math:`h_c(s) = h((s-1)/c)` with :math:`h(z) = z(2-z)` for :math:`0 \leq z < 1` and :math:`h(z) = 1` for :math:`z \geq 1`. These mollifications are :math:`C^1` smooth and ensure that when one contact type vanishes due to mollification, another contact type becomes active, preserving the barrier property.

Usage
-----

GCP is implemented as separate collision and potential classes. A basic example of computing the GCP potential is as follows.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            #include <ipc/ipc.hpp>

            // Set up parameters
            double dhat = 1e-3;     // barrier activation distance
            double alpha_t = 0.5;   // local minimum constraint support
            double beta_t = 0.0;    // local minimum constraint offset
            double alpha_n = 0.1;   // exterior direction constraint support
            double beta_n = 0.0;    // exterior direction constraint offset
            int r = 2;              // barrier exponent (dimension - 1)

            ipc::SmoothContactParameters params(dhat, alpha_t, beta_t, alpha_n, beta_n, r);

            // Build collision set
            bool use_adaptive_dhat = true;
            ipc::SmoothCollisions collisions;
            if (use_adaptive_dhat)
                collisions.compute_adaptive_dhat(collision_mesh, vertices, params);
            collisions.build(collision_mesh, vertices, params, use_adaptive_dhat);

            // Compute potential
            ipc::SmoothContactPotential barrier_potential(params);
            double b = barrier_potential(collisions, collision_mesh, vertices);

            // Compute gradient
            Eigen::VectorXd grad = barrier_potential.gradient(collisions, collision_mesh, vertices);

            // Compute Hessian
            Eigen::SparseMatrix<double> hess = barrier_potential.hessian(
                collisions, collision_mesh, vertices);

    .. md-tab-item:: Python

        .. code-block:: python

            import ipctk

            # Set up parameters
            dhat = 1e-3       # barrier activation distance
            alpha_t = 0.5     # local minimum constraint support
            beta_t = 0.0      # local minimum constraint offset
            alpha_n = 0.1     # exterior direction constraint support
            beta_n = 0.0      # exterior direction constraint offset
            r = 2             # barrier exponent (dimension - 1)

            params = ipctk.SmoothContactParameters(dhat, alpha_t, beta_t, alpha_n, beta_n, r)

            # Build collision set
            use_adaptive_dhat = True
            collisions = ipctk.SmoothCollisions()
            if use_adaptive_dhat:
                collisions.compute_adaptive_dhat(collision_mesh, vertices, params)
            collisions.build(collision_mesh, vertices, params, use_adaptive_dhat)

            # Compute potential
            barrier_potential = ipctk.SmoothContactPotential(params)
            b = barrier_potential(collisions, collision_mesh, vertices)

            # Compute gradient
            grad = barrier_potential.gradient(collisions, collision_mesh, vertices)

            # Compute Hessian
            hess = barrier_potential.hessian(collisions, collision_mesh, vertices)

.. important::
    If ``use_adaptive_dhat`` is true, make sure to call ``SmoothCollisions::compute_adaptive_dhat()`` **before** ``SmoothCollisions::build()``. Adaptive :math:`\hat{d}` computes per-element barrier extents based on the rest configuration to guarantee zero potential (and zero forces) in the undeformed state.

.. note::
    Unlike ``NormalCollisions`` in IPC, ``SmoothCollisions`` must be rebuilt whenever vertex positions change, because the interaction set depends on the current geometry (normals, tangents, and distances).

Parameter Choices
-----------------

The ``SmoothContactParameters`` structure contains the following parameters:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``dhat`` (:math:`\hat{d}`)
     - Global maximum barrier activation distance. The potential vanishes for distances exceeding :math:`\hat{d}`. Unlike IPC, this can be chosen independently of the mesh resolution.
   * - ``alpha_t`` (:math:`\alpha_t`)
     - Controls the support of the smooth Heaviside for the local minimum constraint. Larger values include more points in the interaction set and produce smoother potentials; smaller values are more selective but less smooth.
   * - ``beta_t`` (:math:`\beta_t`)
     - Controls the offset of the smooth Heaviside for the local minimum constraint.
   * - ``alpha_n`` (:math:`\alpha_n`)
     - Controls the support of the smooth Heaviside for the exterior direction constraint.
   * - ``beta_n`` (:math:`\beta_n`)
     - Controls the offset of the smooth Heaviside for the exterior direction constraint.
   * - ``r``
     - The barrier exponent. Should be :math:`n - 1` where :math:`n` is the scene dimension (so :math:`r = 1` for 2D and :math:`r = 2` for 3D).

The :math:`\alpha` and :math:`\beta` parameters for each constraint must satisfy:

.. math::
   0 \leq -\beta < \alpha \leq 1, \quad \alpha + \beta > 0.

As :math:`\alpha + \beta` decreases, the support of the Heaviside function shrinks and the function becomes less smooth, making the potential harder to optimize.

.. tip::
    For simplicity, we recommend setting :math:`\beta_n = \beta_t = 0`. The recommended parameter choices are :math:`\alpha_t \in [0.2, 0.9]` and :math:`\alpha_n = 0.1`. The barrier property is satisfied for any :math:`\alpha > 0`, so the choice only affects smoothness and performance.

Additional internal parameters that may affect behavior:

- **Adaptive dhat ratio** (default ``0.5``): Controls the ratio :math:`\epsilon(x) / d_c(x, f_0)` in the adaptive barrier localization. Can be set via ``SmoothContactParameters::set_adaptive_dhat_ratio()``.
- **Element measure** :math:`L`: For vertices, this is set to the average edge length around the vertex; for edges, to the edge length; for faces, :math:`L` is not needed. This determines the strength of the potential for low-dimensional contact (edge–edge, edge–vertex, vertex–vertex).

Friction
--------

We implement the original friction formulation of :cite:t:`Li2020IPC` following the same style as GCP. Details are covered in the paper :cite:t:`Huang2025GCP`.

Comparison with IPC
-------------------

The following summarizes the key differences between GCP and the original IPC barrier:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Property
     - IPC
     - GCP
   * - Distance criterion
     - Pure Euclidean distance
     - Distance + geometric interaction set
   * - Spurious forces at rest
     - Can occur when :math:`\hat{d}` is large
     - Zero by construction (adaptive :math:`\epsilon`)
   * - Thin shell handling
     - May repel opposite sides
     - Exterior direction constraint eliminates this
   * - :math:`\hat{d}` dependence on mesh
     - Must be small relative to mesh features
     - Independent of discretization
   * - Smoothness
     - :math:`C^2` in distance
     - :math:`C^1` overall (piecewise :math:`C^2`)
   * - Collision rebuild
     - Not needed per step
     - Required each time vertices change

