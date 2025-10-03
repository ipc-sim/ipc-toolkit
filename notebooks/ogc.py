import pathlib
from find_ipctk import ipctk
import polyscope as ps
from polyscope import imgui
import meshio
import numpy as np

import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "mesh",
    type=pathlib.Path,
    help="Path to the mesh file (e.g., .ply)"
)
args = parser.parse_args()

mesh = meshio.read(args.mesh)

# Rotate 90 degrees around the x-axis
R = np.array([[1, 0, 0],
              [0, 0, -1],
              [0, 1, 0]])

cmesh = ipctk.CollisionMesh(
    # mesh.points @ R.T,
    mesh.points,
    ipctk.edges(mesh.cells_dict["triangle"]),
    mesh.cells_dict["triangle"],
)

min_edge_length = np.min(
    np.linalg.norm(cmesh.rest_positions[cmesh.edges[:, 0]] -
                   cmesh.rest_positions[cmesh.edges[:, 1]], axis=1))

max_dhat = 2e-1
dhat = 2e-3  # 0.5 * min_edge_length
use_ogc = True
use_area_weighting = True

# Prebuild the candidates since dhat is the max it will ever be
candidates = ipctk.Candidates()
candidates.build(cmesh, cmesh.rest_positions,
                 max_dhat, broad_phase=ipctk.BVH())

collisions = ipctk.NormalCollisions()


def contact_force():
    global collisions
    collisions.clear()
    collisions.collision_set_type = (
        ipctk.NormalCollisions.CollisionSetType.OGC if use_ogc else ipctk.NormalCollisions.CollisionSetType.IPC)
    collisions.use_area_weighting = use_area_weighting
    collisions.build(candidates, cmesh, cmesh.rest_positions, dhat)

    B = ipctk.BarrierPotential(dhat)
    f = -B.gradient(collisions, cmesh,
                    cmesh.rest_positions).reshape(-1, 3, order="C")
    f *= 1e5

    return f


ps.init()
ps.set_give_focus_on_show(True)
ps.set_ground_plane_mode("shadow_only")

ps_mesh = ps.register_surface_mesh(
    args.mesh.stem,
    cmesh.rest_positions,
    cmesh.faces,
    edge_width=1.0,
)

f = contact_force()
ps_force = ps_mesh.add_scalar_quantity(
    "force",
    np.linalg.norm(f, axis=1),
    defined_on="vertices",
    cmap="reds",
    vminmax=(0, 1e-1),
    enabled=True,
)

ps_force_vector = ps_mesh.add_vector_quantity(
    "force_vector",
    f,
    defined_on="vertices",
    radius=0.01,
    enabled=True,
    # vectortype="ambient",
)

# xi = 0
ev_i = 0
ps_ev = ps.register_curve_network(
    "edge-vertex",
    cmesh.rest_positions[0:3],  # Temporarily use first 3 vertices
    np.array([[1, 2]]),
    enabled=False,
)
ee_i = 0
ps_ee = ps.register_curve_network(
    "edge-edge",
    cmesh.rest_positions[0:4],  # Temporarily use first 4 vertices
    np.array([[0, 1], [2, 3]]),
    enabled=True,
)
ps_ee.set_radius(dhat/2, relative=False)
ps_ee_closest_points = ps.register_point_cloud(
    "edge-edge-closest-points",
    cmesh.rest_positions[0:2],  # Temporarily use first 2 vertices
    enabled=True,
)


def callback():
    global dhat, use_ogc, ev_i, ee_i, max_dhat
    changed0, dhat = imgui.SliderFloat("dhat", dhat, 1e-3, max_dhat, "%.5f")
    if changed0 and dhat > max_dhat:
        max_dhat = dhat
        candidates.build(cmesh, cmesh.rest_positions, max_dhat,
                         broad_phase=ipctk.BVH())
    changed1, use_ogc = imgui.Checkbox("use_ogc", use_ogc)
    if changed0 or changed1:
        f = contact_force()
        ps_mesh.add_scalar_quantity("force", np.linalg.norm(f, axis=1))
        ps_mesh.add_vector_quantity("force_vector", f)
    imgui.Text(
        f"Number of collisions: {len(collisions.vv_collisions)} VV, {len(collisions.ev_collisions)} EV,"
        f" {len(collisions.ee_collisions)} EE, {len(collisions.fv_collisions)} FV")

    if len(collisions.ev_collisions) > 0:
        ev_i = min(ev_i, len(collisions.ev_collisions) - 1)
        _, ev_i = imgui.SliderInt("Edge-Vertex Index", ev_i,
                                  0, len(collisions.ev_collisions) - 1)
        ps_ev.update_node_positions(cmesh.rest_positions[
            collisions.ev_collisions[ev_i].vertex_ids(cmesh.edges, cmesh.faces)[:3]])
    else:
        ps_ev.set_enabled(False)

    if len(collisions.ee_collisions) > 0:
        ee_i = min(ee_i, len(collisions.ee_collisions) - 1)
        changed2, ee_i = imgui.SliderInt("Edge-Edge Index", ee_i,
                                         0, len(collisions.ee_collisions) - 1)
        coeffs = collisions.ee_collisions[ee_i].compute_coefficients(
            cmesh.rest_positions, cmesh.edges, cmesh.faces)
        imgui.Text("Coeffs: {}".format(coeffs))
        if changed0 or changed1 or changed2:
            V_ee = cmesh.rest_positions[
                collisions.ee_collisions[ee_i].vertex_ids(cmesh.edges, cmesh.faces)]
            ps_ee.update_node_positions(V_ee)
            if changed0:
                ps_ee.set_radius(dhat/2, relative=False)
            ps_ee_closest_points.update_point_positions(
                np.vstack([
                    coeffs[0] * V_ee[0] + coeffs[1] * V_ee[1],
                    -coeffs[2] * V_ee[2] - coeffs[3] * V_ee[3]]))
    else:
        ps_ee.set_enabled(False)
        ps_ee_closest_points.set_enabled(False)

    # if changed2:
    #     adj_vi = cmesh.edge_vertex_adjacencies[ei]
    #     if len(adj_vi) == 1:
    #         adj_vi = [adj_vi[0], adj_vi[0]]
    #     ps_edge.update_node_positions(
    #         cmesh.rest_positions[[
    #             *cmesh.edges[ei],
    #             *adj_vi
    #         ]])
    #
    # _, xi = imgui.SliderInt("Vertex Index##mine", xi, 0,
    #                         len(cmesh.rest_positions) - 1)
    # x = cmesh.rest_positions[xi]
    # imgui.Text(f"Vertex position: {x}")
    # x0, x1 = cmesh.rest_positions[cmesh.edges[ei]]
    # passed = True
    # for vi in cmesh.edge_vertex_adjacencies[ei]:
    #     x2 = cmesh.rest_positions[vi]
    #     alpha = ipctk.point_edge_closest_point(x2, x0, x1)
    #     imgui.Text(f"alpha: {alpha:.5f}")
    #     p = (x1 - x0) * alpha + x0
    #     imgui.Text(f"p: {p}, x2: {x2}, x: {x}")
    #     imgui.Text(f"b: {(x - p).dot(p - x2)}")
    #     passed = passed and (x - p).dot(p - x2) > 0
    # imgui.Text(f"check_edge_feasible_region: {passed}")


ps.set_user_callback(callback)


ps.show()
