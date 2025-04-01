import numpy as np
import scipy as sp

import find_ipctk
from ipctk import CollisionMesh, NormalCollisions, TangentialCollisions, FrictionPotential, BarrierPotential

def test_collision_mesh_material_ids():
    # Create a simple mesh
    vertices = np.array([
        [0, 0, 0],
        [1, 0, 0], 
        [0, 1, 0],
        [1, 1, 0]
    ], dtype=float)
    
    edges = np.array([
        [0, 1],
        [1, 3],
        [3, 2],
        [2, 0]
    ], dtype=int)
    
    # Create a collision mesh
    mesh = CollisionMesh(vertices, edges)
    
    # Set material IDs for vertices
    vertex_materials = np.array([0, 0, 1, 1], dtype=int)
    mesh.vertex_materials = vertex_materials
    
    # Check material IDs were set correctly
    assert mesh.has_material_ids()
    assert np.all(mesh.vertex_materials == vertex_materials)
    
    # Set a single material ID for everything
    mesh.set_single_material_id(2)
    assert np.all(mesh.vertex_materials == 2)
    assert np.all(mesh.edge_materials == 2)
    
    # Test NO_MATERIAL_ID case
    # Create a mesh without setting materials
    mesh2 = CollisionMesh(vertices, edges)
    assert not mesh2.has_material_ids()

def test_material_friction():
    # Create a simple scene with two objects of different materials
    vertices = np.array([
        [0, 0, 0],
        [1, 0, 0], 
        [0, 1, 0],
        [1, 1, 0]
    ], dtype=float)
    
    edges = np.array([
        [0, 1],
        [1, 3],
        [3, 2],
        [2, 0]
    ], dtype=int)
    
    # Create a collision mesh
    mesh = CollisionMesh(vertices, edges)
    
    # Set different materials for vertices
    vertex_materials = np.array([0, 0, 1, 1], dtype=int)
    mesh.vertex_materials = vertex_materials
    
    # Create a material-specific friction table
    friction_table = np.array([
        [0.3, 0.5],  # material 0 against 0 or 1
        [0.5, 0.7]   # material 1 against 0 or 1
    ])
    
    # Create normal collisions
    normal_collisions = NormalCollisions()
    normal_collisions.build(mesh, vertices, 0.1)
    
    # Create tangential collisions with material-specific friction
    tangential_collisions = TangentialCollisions()
    tangential_collisions.build(
        mesh, vertices, normal_collisions, 
        BarrierPotential(0.1), 100.0, 0.3  # default mu is 0.3
    )
    
    # Create a friction potential and set the material friction table
    friction_potential = FrictionPotential(1e-5)
    friction_potential.set_material_friction_table(friction_table)
    
    # Apply velocities to test friction
    velocities = np.zeros_like(vertices)
    velocities[2:4, 0] = 0.1  # Move vertices with material 1
    
    # Compute friction forces
    forces = friction_potential.force(
        tangential_collisions, mesh, vertices, 
        np.zeros_like(vertices), velocities,
        BarrierPotential(0.1), 100.0
    )
    
    # The forces should be non-zero where there are collisions
    assert np.linalg.norm(forces) > 0
