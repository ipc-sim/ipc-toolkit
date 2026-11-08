Rigid and Affine Body Dynamics
==============================

The toolkit can simulate bodies that move as a single unit instead of deforming like a finite element mesh. It supports two models:

* **Rigid bodies**, following Rigid IPC :cite:t:`Ferguson2021RigidIPC`. Each body has 6 degrees of freedom (DOFs) in 3D: a position and a rotation vector. Bodies never deform, and CCD follows their curved trajectories (see :doc:`nonlinear_ccd`).
* **Affine bodies**, following Affine Body Dynamics (ABD) `[Lan et al. 2022] <https://doi.org/10.1145/3528223.3530064>`__. Each body has 12 DOFs in 3D: a position and a :math:`3 \times 3` matrix :math:`A`. A stiff penalty keeps :math:`A` close to a rotation, and each Newton step moves the vertices along straight lines, so we can use ordinary CCD. Affine bodies also support :ref:`joints <dynamics-joints>`.

Both models run the same code. We write the body terms (inertia, gravity and other forces, orthogonality, and so on) once in affine coordinates :math:`y = [p; \operatorname{vec}(A)]`. For rigid bodies, :cpp:class:`ipc::dynamics::RigidToAffine` maps the DOFs :math:`[p; \theta]` to :math:`A = R(\theta)` and applies the chain rule, while :cpp:class:`ipc::dynamics::AffineToAffine` is the identity.

.. note::
    The library only provides the math: potentials, the change of variables, time integration (:cpp:class:`ipc::dynamics::BDF`), and CCD. It does not include a solver. The **demo simulator** (``ipc_toolkit_demo``) combines these parts with `PolySolve <https://github.com/polyfem/polysolve>`_'s Newton solver into a complete stepping loop. This tutorial uses it. Everything in the demo uses public library components, so you can copy it as a starting point for your own simulator.

Dropping a Cube
---------------

The demo is built by default for top-level builds (CMake option ``IPC_TOOLKIT_BUILD_DEMO``), and the Python bindings expose it as ``ipctk.demo``. Given a triangle mesh with vertices ``V``, edges ``E``, and faces ``F``, the following code drops a cube onto the ground:

.. md-tab-set::

    .. md-tab-item:: C++

        .. literalinclude:: ../../../tests/src/tests/demo/test_tutorial.cpp
            :language: c++
            :start-after: // BEGIN_DROP_A_CUBE
            :end-before: // END_DROP_A_CUBE
            :dedent: 4

    .. md-tab-item:: Python

        .. literalinclude:: ../../../python/tests/test_dynamics.py
            :language: python
            :start-after: # BEGIN_DROP_A_CUBE
            :end-before: # END_DROP_A_CUBE
            :dedent: 8

Let's walk through it:

#. :cpp:func:`ipc::rigid::RigidBodies::build_from_meshes` takes one mesh per body, along with a density and an initial pose for each. It re-centers each body on its center of mass and computes its mass properties. All the bodies live in one :cpp:class:`ipc::CollisionMesh`, so the toolkit's collision, friction, and adhesion potentials work on them unchanged.
#. ``planes`` holds infinite planes. They collide like any other surface but never move.
#. ``Simulator::Settings`` collects the simulation parameters. The defaults give a good first result: gravity is :math:`(0, -9.81, 0)`, the barrier activation distance ``dhat`` is 1 mm, and the time integrator is implicit Euler (``bdf_order = 1``).
#. ``Simulator::run`` steps the simulation to ``t_end``. Use ``step`` to advance one step at a time instead.
#. ``rigid_poses`` returns the current pose of every body. Pass the poses to ``RigidBodies::vertices`` to get world-space vertices, or use ``rigid_pose_history`` for every step.

Choosing the Model
------------------

The only change needed to simulate affine bodies is ``settings.body_dynamics``, which the code above already sets (``RIGID`` or ``AFFINE``). Here is how they differ:

.. list-table::
    :header-rows: 1

    * -
      - Rigid
      - Affine
    * - DOFs per body (3D / 2D)
      - 6 / 3
      - 12 / 6
    * - CCD
      - Nonlinear, along curved trajectories
      - Linear
    * - Shape
      - Exactly rigid
      - Nearly rigid; stretching is possible but penalized by ``orthogonality_stiffness``
    * - Joints
      - No
      - Yes

Friction and Kinematic Bodies
-----------------------------

Set ``settings.friction_coefficient`` above zero to enable lagged friction :cite:p:`Ferguson2021RigidIPC`. ``friction_iterations`` controls how many times the friction is re-lagged per step.

Each body also has a type, ``RigidBody::Type``, which defaults to ``DYNAMIC``. A ``STATIC`` body never moves but still collides. A ``KINEMATIC`` body follows a prescribed motion, and pushes dynamic bodies out of the way. Set the type before you build the simulator:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            (*bodies)[0].set_type(rigid::RigidBody::Type::STATIC);

    .. md-tab-item:: Python

        .. code-block:: python

            bodies[0].type = ipctk.RigidBody.Type.STATIC

A kinematic body moves with its initial velocity by default. To script its motion, attach a ``KinematicDriver`` with ``Simulator::set_kinematic_driver``.

.. _dynamics-joints:

Joints
------

Affine bodies can be connected with linear equality constraints on their DOFs, as in `Chen et al. 2022, §4.1 <https://doi.org/10.1145/3528223.3530076>`__. :cpp:class:`ipc::affine::JointConstraints` provides ``add_point_connection``, ``add_hinge``, ``add_sliding_plane``, ``add_fixed_point``, and ``add_fixed_body``. Pass the constraints to the ``Simulator`` constructor, and it enforces them exactly by a change of variables to a reduced set of coordinates. See ``python/examples/joints.py`` for a double pendulum.

Reading and Writing Scenes
--------------------------

When the toolkit is built with ``IPC_TOOLKIT_WITH_GLTF=ON``, you can load a scene from a glTF file and save a simulation as an animation:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            #include <ipc/io/read_gltf.hpp>
            #include <ipc/io/write_gltf.hpp>

            auto [bodies, poses] = rigid::read_gltf("scene.glb", /*convert_planes=*/true);
            // ... build a Simulator and run it ...
            rigid::write_gltf("animation.glb", *bodies, sim.rigid_pose_history(), dt);

    .. md-tab-item:: Python

        .. code-block:: python

            bodies, poses = ipctk.read_gltf("scene.glb", convert_planes=True)
            # ... build a Simulator and run it ...
            ipctk.write_gltf("animation.glb", bodies, sim.rigid_pose_history, dt)

Each mesh node becomes a body. ``convert_planes`` turns any mesh made of two coplanar triangles into an infinite plane, which is a convenient way to author the ground. The file only holds geometry: set densities, body types, and initial velocities in code. ``python/examples/rigid.py`` is an interactive viewer that does all of this.

.. tip::
    PolySolve logs through its own logger, separate from the toolkit's. Use ``Settings.solver_log_level`` to control its verbosity.

Going Further
-------------

* The `demo README <https://github.com/ipc-sim/ipc-toolkit/tree/main/demo>`__ explains the design of the simulator in detail: the incremental potential and its scaling, how rigid and affine bodies share one formulation, and how the solver differs from the original Rigid IPC.
* To build your own stepping loop, combine :cpp:class:`ipc::dynamics::BodyPotentials` (sums the body terms and applies the chain rule once), :cpp:class:`ipc::dynamics::BDF` (time integration), the barrier and friction potentials, and a solver of your choice. See the :doc:`../cpp-api/dynamics` and :doc:`../cpp-api/rigid_bodies` API pages.
