Rigid Bodies
============

Rigid body dynamics, following Rigid IPC :cite:`Ferguson2021RigidIPC`. Each body
is re-centered on its center of mass and expressed in its principal frame, so
the mass matrix is diagonal. The bodies of a scene are stored together as one
:cpp:class:`ipc::CollisionMesh`, which lets every contact, friction, and
adhesion potential work on rigid bodies unchanged.

Rigid Bodies
------------

.. doxygenclass:: ipc::rigid::RigidBodies
    :allow-dot-graphs:

Rigid Body
----------

.. doxygenclass:: ipc::rigid::RigidBody
    :allow-dot-graphs:

Pose
----

.. doxygenstruct:: ipc::rigid::Pose
    :allow-dot-graphs:

Rotation Parameterization
^^^^^^^^^^^^^^^^^^^^^^^^^

A rigid body's rotation is a rotation vector in 3D and an angle in 2D. These
functions convert it to a rotation matrix and provide the derivatives needed by
the chain rule.

.. doxygenfunction:: ipc::rigid::rotation_angle_to_matrix
.. doxygenfunction:: ipc::rigid::rotation_vector_to_matrix
.. doxygenfunction:: ipc::rigid::rotation_vector_to_matrix_jacobian
.. doxygenfunction:: ipc::rigid::rotation_vector_to_matrix_hessian
.. doxygenfunction:: ipc::rigid::rotation_to_matrix_jacobian
.. doxygenfunction:: ipc::rigid::rotation_to_matrix_hessian
.. doxygenfunction:: ipc::rigid::rotation_matrix_to_vector

Mass Properties
---------------

.. doxygenfunction:: ipc::rigid::compute_mass_properties
.. doxygenfunction:: ipc::rigid::construct_mass_matrix

Rigid Trajectory
----------------

The trajectory of a point on a rigid body is curved, not a line. Nonlinear CCD
(see :doc:`ccd`) follows it piece by piece, using the chord deviation to bound
how far each piece strays from its straight-line approximation.

.. doxygenclass:: ipc::rigid::RigidTrajectory
    :allow-dot-graphs:

.. doxygenfunction:: ipc::rigid::rigid_chord_deviation
.. doxygenfunction:: ipc::rigid::relative_rigid_chord_deviation

Rigid Candidates
----------------

.. doxygenclass:: ipc::rigid::RigidCandidates
    :allow-dot-graphs:

Contact Hessian
---------------

.. doxygenfunction:: ipc::rigid::per_contact_rigid_hessian
