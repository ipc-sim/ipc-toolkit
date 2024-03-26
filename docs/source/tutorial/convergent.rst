.. _convergent-collision-formulation:

Convergent Formulation
======================

In addition to the original implementation of :cite:t:`Li2020IPC`, we also implement the convergent formulation of :cite:t:`Li2023Convergent`.

To enable the convergent formulation, we need to set ``use_convergent_formulation`` before building ``Collisions``:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            collisions.set_use_convergent_formulation(true);
            collisions.build(collision_mesh, vertices, dhat);

    .. md-tab-item:: Python

        .. code-block:: python

            collisions.use_convergent_formulation = True
            collisions.build(collision_mesh, vertices, dhat)

.. important::
    The variable ``use_convergent_formulation`` should be set before calling ``Collisions::build`` for it to take effect. By default, it is ``false``.

Technical Details
-----------------

*We briefly summarize the convergent formulation here for convenience.*

In order to derive a convergent formulation, we first define a continuous form of our barrier potential :math:`P`. For a surface :math:`\mathcal{S}` embedded in 3D space, we parameterize the surfaces by common (possibly discontinuous) coordinates :math:`u \in \tilde{M} \subset \mathbb{R}^2`, so that :math:`\mathcal{S}(u)` traverses the material points across all surfaces contiguously. The total barrier potential is then

.. math::
   P(\mathcal{S})=\frac{1}{2} \int_{u \in \tilde{M}} \max _{v \in \tilde{M} \setminus{ }_r u} b(d(\mathcal{S}(u), \mathcal{S}(v)), \hat{d})~\mathrm{d} u,

where we define the operator :math:`\setminus_r: \mathcal{P}(\mathbb{R}^2) \times \mathbb{R} \times \mathbb{R}^2 \mapsto \mathcal{P}(\mathbb{R}^2)` to be

.. math::
    \tilde{M} \setminus_r u:=\left\{v \in \tilde{M} \mid\|u-v\|_2>r\right\}

with :math:`r \rightarrow 0`.

We then define our surface discretization with a triangulated boundary mesh geometry. As in the smooth case, we can parameterize the domain across all triangles with :math:`u \in \tilde{M}` so that :math:`p(u): \tilde{M} \mapsto \mathbb{R}^3` traverses all material points, across all triangles :math:`t ∈ T` in the triangle mesh contiguously. The corresponding surface barrier potential is then

.. math::
    \frac{1}{2} \int_{u \in \tilde{M}} \max_{t \in T \backslash p(u)} b(d(p(u), t), \hat{d})~\mathrm{d} u,

where :math:`T \setminus p` is the set of boundary faces that do not contain the point :math:`p`.

Applying mesh vertices as nodes (and quadrature points), we numerically integrate the surface barrier potential. For each nodal position :math:`x \in V` we then have a corresponding material space coordinate :math:`\bar{x} \in \bar{V}`. Piecewise linear integration of the surface barrier is then

.. math::
    \frac{1}{2} \sum_{\bar{x} \in \bar{V}} w_{\bar{x}} \max _{t \in T \backslash x(\bar{x})} b(d(x(\bar{x}), t), \hat{d}),

where :math:`w_{\bar{x}}` are the quadrature weights, each given by one-third of the sum of the areas (in material space) of the boundary triangles incident to :math:`\bar{x}`.

We next need to smoothly approximate the max operator in the barrier potentials. However, common approaches such as an :math:`L^p`-norm or LogSumExp would decrease sparsity in subsequent numerical solves by increasing the stencil size per collision evaluation. We instead leverage the locality of our barrier function to approximate the max operator by removing duplicate distance pairs. Our resulting approximators for a triangulated surface is

.. math::
    \begin{aligned}
    \Psi_s(x) & =\sum_{t \in T \backslash x} b(d(x, t), \hat{d})-\sum_{e \in E_{\text {int }} \backslash x} b(d(x, e), \hat{d})+\sum_{x_2 \in V_{i n t} \backslash x} b\left(d\left(x, x_2\right), \hat{d}\right) \\
    & \approx \max _{t \in T \backslash x} b(d(x, t), \hat{d}),
    \end{aligned}

where :math:`V_{\text{int}} \subseteq V` is the subset of internal surface nodes and :math:`E_{\text{int}} \subseteq E` is the subset of internal surface edges (i.e., edges incident to two triangles). For locally convex regions this estimator is tight while remaining smooth. In turn, for nonconvex regions, it improves over direct summation.

The corresponding discrete barrier potential is then simply

.. math::
    P_s(V)= \frac{1}{2} \sum_{x \in V} w_x \Psi_s(x),

where we simplify with :math:`w_x = w_{\bar{x}}` defined appropriately, per domain, as covered above.

Please, see the `paper <https://arxiv.org/abs/2307.15908>`_ for more details (including the formulation for 2D curves and edge-edge collisions) and evaluation.

The key difference between the original and the convergent formulations is that we (1) include area weights in the barrier potential and (2) include additional (negative) terms to cancel out the duplicate distance pairs. While this requires minor algorithmic changes, it produces considerably better results.

Physical Barrier
----------------

We want our barrier potential to have the same units as our elastic potential (e.g., :math:`\text{J}`). Together with the area weighting (discussed above), this means the barrier should have units of pressure times distance (e.g., :math:`\text{Pa} \cdot \text{m}`). That is,

.. math::
    \text{Pa} \cdot \text{m} \cdot \text{m}^2 = \frac{\text{N}}{\text{m}^2} \cdot \text{m} \cdot \text{m}^2 = \text{N} \cdot \text{m} = \text{J}.

To achieve this, (when using the convergent formulation) we modify the barrier function to have units of distance:

.. math::
    b(d, \hat{d})=\left\{\begin{array}{lr}
    -\hat{d}\left(\frac{d}{\hat{d}}-1\right)^2 \ln \left(\frac{d}{\hat{d}}\right), & 0<d<\hat{d} \\
    0 & d \geq \hat{d}
    \end{array}\right.

.. note::
    This is equivalent to the original barrier function of :cite:p:`Li2020IPC` times :math:`1/\hat{d}^3` when using squared distances. Therefore, to simplify the implementation we only implement the original barrier function and multiply all barrier potentials by :math:`1/\hat{d}^3`.

The barrier stiffness (:math:`\kappa`) then has units of pressure (e.g., :math:`\text{Pa}`), the same as Young's modulus (:math:`E`) in elasticity.
This implies we can get good solver convergence even when using a fixed :math:`\kappa` by setting it relative to the material's Young's modulus (:math:`\kappa = 0.1 E` works well in many examples).
The intention is to treat the barrier as a thin elastic region around the mesh, and having consistent units makes it easier to pick the stiffness for this "material".

.. _convergent-friction-formulation:

Friction
--------

Just as with the :ref:`collisions <convergent-collision-formulation>`, we implement both the original friction formulation of :cite:t:`Li2020IPC` and the convergent formulation of :cite:t:`Li2023Convergent`.

The choice of formulation is dependent on how the fixed set of ``collisions`` given to ``FrictionCollisions::build`` was built. If the ``collisions`` were built using the convergent formulation, then the friction collisions will also use the convergent formulation. Otherwise, the original formulation will be used.
