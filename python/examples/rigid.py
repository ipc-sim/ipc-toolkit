from find_ipctk import ipctk
import meshio
import polyscope as ps
from polyscope import imgui
import numpy as np

import pathlib

mesh = meshio.read(pathlib.Path(__file__).parents[2] / "tests/data/cube.ply")

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

for pose in initial_poses:
    print(pose)

ps.init()

ps.set_give_focus_on_show(True)

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

playing = False


def foo():
    global playing
    if imgui.Button("Play" if not playing else "Pause"):
        playing = not playing
    imgui.SameLine()
    if imgui.Button("Step") or playing:
        sim.step()
        ps_mesh.update_vertex_positions(
            bodies.vertices(sim.poses_history[-1])
        )
    imgui.SameLine()
    if imgui.Button("Reset"):
        sim.reset()
        ps_mesh.update_vertex_positions(
            bodies.vertices(sim.poses_history[-1])
        )


ps.set_user_callback(foo)

ps.show()
