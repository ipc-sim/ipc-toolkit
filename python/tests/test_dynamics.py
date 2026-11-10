import numpy as np
from find_ipctk import ipctk
from utils import load_mesh


def test_dropped_cube_comes_to_rest_on_the_ground():
    if not hasattr(ipctk, "demo"):
        return  # The demo simulator is not part of this build.

    # A unit cube centered at the origin with vertices V, edges E, and faces F.
    V, E, F = load_mesh("cube.ply")

    for body_dynamics in (
        ipctk.demo.Simulator.BodyDynamics.RIGID,
        ipctk.demo.Simulator.BodyDynamics.AFFINE,
    ):
        # BEGIN_DROP_A_CUBE
        # Hold the cube 0.5 m above the ground (its center is at y = 1).
        poses = ipctk.Poses([ipctk.Pose(np.array([0.0, 1.0, 0.0]), np.zeros(3))])

        bodies = ipctk.RigidBodies(
            rest_positions=[V],
            edges=[E],
            faces=[F],
            densities=[1000.0],
            initial_poses=poses,
        )

        # The ground is the plane y = 0. `planes` returns a copy, so we assign
        # the whole list.
        bodies.planes = [ipctk.Hyperplane(np.array([0.0, 1.0, 0.0]), np.zeros(3))]

        settings = ipctk.demo.Simulator.Settings()
        settings.body_dynamics = body_dynamics  # RIGID or AFFINE

        sim = ipctk.demo.Simulator(bodies, poses, dt=0.01, settings=settings)
        success = sim.run(0.6)

        final_poses = sim.rigid_poses
        # END_DROP_A_CUBE

        assert success

        # The cube has landed: it is on the ground (not through it), and it has
        # fallen at least most of the way down.
        height = final_poses[0].position[1]
        assert height >= 0.5 - 1e-8
        assert height < 0.6
