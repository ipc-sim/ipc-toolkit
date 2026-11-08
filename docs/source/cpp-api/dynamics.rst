Dynamics
========

Components for assembling a body dynamics simulation. The library provides the
math (potentials, change of variables, time integration) and leaves the solver
and the stepping loop to the caller. See the `demo simulator
<https://github.com/ipc-sim/ipc-toolkit/tree/main/demo>`__ for an example.

Change of Variables
-------------------

A :cpp:class:`ipc::dynamics::ToAffine` maps a body's degrees of freedom to its
affine pose, so the same body terms (see :doc:`affine_bodies`) serve both
rigid and affine bodies.

.. doxygenclass:: ipc::dynamics::ToAffine
    :allow-dot-graphs:

.. doxygenclass:: ipc::dynamics::RigidToAffine
    :allow-dot-graphs:

.. doxygenclass:: ipc::dynamics::AffineToAffine
    :allow-dot-graphs:

Body Potentials
---------------

.. doxygenclass:: ipc::dynamics::BodyPotentials
    :allow-dot-graphs:

Time Integration
----------------

.. doxygenclass:: ipc::dynamics::ImplicitTimeIntegrator
    :allow-dot-graphs:

.. doxygenclass:: ipc::dynamics::BDF
    :allow-dot-graphs:
