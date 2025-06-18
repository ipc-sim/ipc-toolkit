from find_ipctk import ipctk
import meshio
import polyscope as ps
import numpy as np

import pathlib

mesh = meshio.read(pathlib.Path(__file__).parents[2] / "tests/data/bunny.ply")

initial_poses = ipctk.Poses([
    ipctk.Pose(position=np.zeros(3), rotation=np.zeros(3))
])

bodies = ipctk.RigidBodies(
    rest_positions=[mesh.points.astype("float64")],
    edges=[ipctk.edges(mesh.cells_dict["triangle"])],
    faces=[mesh.cells_dict["triangle"]],
    densities=[1.0],
    initial_poses=initial_poses
)

ps.init()

ps_mesh = ps.register_surface_mesh(
    "bunny",
    bodies.vertices(initial_poses),
    bodies.faces
)

sim = ipctk.Simulator(
    bodies=bodies,
    initial_poses=initial_poses,
    dt=0.01
)


def foo():
    sim.step()
    ps_mesh.update_vertex_positions(
        bodies.vertices(sim.poses_history[-1])
    )


ps.set_user_callback(foo)

ps.show()
