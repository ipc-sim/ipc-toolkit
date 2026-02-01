import pathlib, meshio
import numpy as np
import polyscope as ps
from find_ipctk import ipctk
from polyscope import imgui
from scipy.spatial.transform import Rotation

mesh = meshio.read(pathlib.Path(__file__).parents[2] / "tests/data/bunny-lowpoly.obj")

initial_poses = ipctk.Poses([ipctk.Pose(position=np.zeros(3), rotation=np.zeros(3))])
initial_poses[0].position = np.array([0, 0.5, 0])  # Lift the bunny above the ground
initial_poses[0].rotation = np.array([np.pi / 4, 0, 0])  # Rotate around x-axis

bodies = ipctk.RigidBodies(
    rest_positions=[mesh.points.astype("float64")],
    edges=[ipctk.edges(mesh.cells_dict["triangle"])],
    faces=[mesh.cells_dict["triangle"]],
    densities=[1.0],
    initial_poses=initial_poses,
)

# for pose in initial_poses:
#     print(pose)

ps.init()

ps.set_ground_plane_height_mode("manual")

ps.set_give_focus_on_show(True)

ps_mesh = ps.register_surface_mesh(
    "rigid body",
    bodies.vertices(initial_poses),
    bodies.faces,
)

ps_com = ps.register_point_cloud(
    "rigid body com",
    initial_poses[0].position.reshape(-1, 3),
)

for i in range(3):
    dim = np.zeros(3)
    dim[i] = 100 * bodies[0].moment_of_inertia[i]
    R = Rotation.from_rotvec(initial_poses[0].rotation.copy()).as_matrix()
    ps_com.add_vector_quantity(
        "xyz"[i],
        (R @ dim).reshape(-1, 3),
        enabled=True,
        vectortype="ambient",
    )

sim = ipctk.Simulator(
    bodies=bodies,
    initial_poses=initial_poses,
    dt=1 / 60.0,
)

playing = False


def callback():
    global playing
    if imgui.Button("Play" if not playing else "Pause") or imgui.IsKeyPressed(
        imgui.ImGuiKey_Space
    ):
        playing = not playing
    imgui.SameLine()
    if imgui.Button("Step") or playing:
        sim.step()
        ps_mesh.update_vertex_positions(bodies.vertices(sim.pose_history[-1]))
    imgui.SameLine()
    if imgui.Button("Reset"):
        sim.reset()
        ps_mesh.update_vertex_positions(bodies.vertices(sim.pose_history[-1]))


ps.set_user_callback(callback)

ps.show()
