Affine Bodies
=============

The body (non-contact) terms of the simulation, written once in affine
coordinates :math:`y = [p; \operatorname{vec}(A)]`. They follow Affine Body
Dynamics (ABD) `[Lan et al. 2022] <https://doi.org/10.1145/3528223.3530064>`__.
Rigid bodies reuse the same terms through the change of variables in
:doc:`dynamics`.

Every term except inertia is a pure, unscaled physical potential (energy
units). The caller applies the time-integration scaling when it sums them.

Affine Pose
-----------

.. doxygenstruct:: ipc::affine::Pose
    :allow-dot-graphs:

.. doxygenfunction:: ipc::affine::dof_to_pose
.. doxygenfunction:: ipc::affine::dof_to_poses
.. doxygenfunction:: ipc::affine::vertices
.. doxygenfunction:: ipc::affine::affine_jacobian

Inertial Term
-------------

.. doxygenclass:: ipc::affine::InertialTerm
    :allow-dot-graphs:

Body Forces
-----------

.. doxygenclass:: ipc::affine::BodyForces
    :allow-dot-graphs:

Orthogonality Potential
-----------------------

A stiff penalty that keeps the matrix :math:`A` close to :math:`SO(\text{dim})`.
It applies to affine bodies only.

.. doxygenclass:: ipc::affine::OrthogonalityPotential
    :allow-dot-graphs:

Augmented Lagrangian
--------------------

Drives kinematic bodies to their target poses.

.. doxygenclass:: ipc::affine::AugmentedLagrangian
    :allow-dot-graphs:

Joint Constraints
-----------------

Linear equality constraints between bodies, as in `Chen et al. 2022, §4.1
<https://doi.org/10.1145/3528223.3530076>`__. They are enforced exactly by a
change of variables to a reduced set of coordinates.

.. doxygenclass:: ipc::affine::JointConstraints
    :allow-dot-graphs:
