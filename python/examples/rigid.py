from find_ipctk import ipctk
import meshio
import polyscope as ps
from polyscope import imgui
import numpy as np
from scipy.spatial.transform import Rotation

import pathlib

mesh = meshio.read(pathlib.Path(
    __file__).parents[2] / "tests/data/wing-nut.obj")

initial_poses = ipctk.Poses([
    ipctk.Pose(position=np.zeros(3), rotation=np.zeros(3))
])
# initial_poses[0].rotation = np.array([np.pi, 0, 0])  # Rotate around x-axis

bodies = ipctk.RigidBodies(
    rest_positions=[mesh.points.astype("float64")],
    edges=[ipctk.edges(mesh.cells_dict["triangle"])],
    faces=[mesh.cells_dict["triangle"]],
    densities=[1.0],
    initial_poses=initial_poses
)

# for pose in initial_poses:
#     print(pose)

ps.init()

ps.set_give_focus_on_show(True)

ps_mesh = ps.register_surface_mesh(
    "rigid body",
    bodies.vertices(initial_poses),
    bodies.faces
)

ps_com = ps.register_point_cloud(
    "rigid body com",
    initial_poses[0].position.reshape(-1, 3)
)

for i in range(3):
    dim = np.zeros(3)
    dim[i] = 100 * bodies[0].moment_of_inertia[i]
    R = Rotation.from_rotvec(
        initial_poses[0].rotation.copy()).as_matrix()
    ps_com.add_vector_quantity(
        "xyz"[i],
        (R @ dim).reshape(-1, 3),
        enabled=True,
        vectortype="ambient"
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
            bodies.vertices(sim.pose_history[-1])
        )
    imgui.SameLine()
    if imgui.Button("Reset"):
        sim.reset()
        ps_mesh.update_vertex_positions(
            bodies.vertices(sim.pose_history[-1])
        )


ps.set_user_callback(foo)

ps.show()
