"""
Anisotropic Friction Example

This example demonstrates how to use anisotropic friction in IPC Toolkit.
Anisotropic friction allows different friction coefficients along different
tangent directions, creating direction-dependent friction behavior.

The example shows:
1. Setting up a simple collision scenario
2. Building normal and tangential collisions
3. Assigning anisotropic friction coefficients
4. Computing effective friction for different velocity directions
"""

import numpy as np
from find_ipctk import ipctk


def create_simple_scene():
    """Create a simple scene with two vertices close together."""
    # Two vertices in 3D space, close enough to create a collision
    vertices = np.array([
        [0.0, 0.0, 0.0],   # Vertex 0
        [0.01, 0.0, 0.0],  # Vertex 1 (very close to vertex 0)
    ])
    
    # No edges or faces needed for vertex-vertex collision
    edges = np.array([], dtype=np.int32).reshape(0, 2)
    faces = np.array([], dtype=np.int32).reshape(0, 3)
    
    return vertices, edges, faces


def main():
    print("=" * 70)
    print("Anisotropic Friction Example")
    print("=" * 70)
    print()
    
    # -------------------------------------------------------------------------
    # Step 1: Setup - Create simple scene
    # -------------------------------------------------------------------------
    print("Step 1: Creating simple scene with two vertices...")
    vertices, edges, faces = create_simple_scene()
    print(f"  Created {len(vertices)} vertices")
    print()
    
    # -------------------------------------------------------------------------
    # Step 2: Build Collisions
    # -------------------------------------------------------------------------
    print("Step 2: Building collision mesh and collisions...")
    
    # Contact parameters
    dhat = 0.01  # Distance threshold for collisions
    barrier_stiffness = 1.0
    mu_s = 0.5  # Initial isotropic static friction
    mu_k = 0.3  # Initial isotropic kinetic friction
    
    # Build collision mesh
    collision_mesh = ipctk.CollisionMesh(vertices, edges, faces)
    print(f"  Collision mesh created with {collision_mesh.num_vertices()} vertices")
    
    # Build normal collisions
    normal_collisions = ipctk.NormalCollisions()
    normal_collisions.build(collision_mesh, vertices, dhat)
    print(f"  Found {normal_collisions.size()} normal collision(s)")
    
    if normal_collisions.size() == 0:
        print("  WARNING: No normal collisions found. Cannot demonstrate friction.")
        return
    
    # Build tangential collisions with initial isotropic friction
    tangential_collisions = ipctk.TangentialCollisions()
    barrier_potential = ipctk.BarrierPotential(dhat)
    tangential_collisions.build(
        collision_mesh, vertices, normal_collisions,
        barrier_potential, barrier_stiffness, mu_s, mu_k
    )
    print(f"  Found {tangential_collisions.size()} tangential collision(s)")
    print()
    
    if tangential_collisions.size() == 0:
        print("  WARNING: No tangential collisions found. Cannot demonstrate friction.")
        return
    
    # -------------------------------------------------------------------------
    # Step 3: Assign Anisotropic Friction Coefficients
    # -------------------------------------------------------------------------
    print("Step 3: Assigning anisotropic friction coefficients...")
    
    # Define anisotropic friction coefficients
    # Higher friction in first tangent direction, lower in second
    mu_s_aniso = np.array([0.8, 0.4])  # Static: 0.8 along t0, 0.4 along t1
    mu_k_aniso = np.array([0.6, 0.3])  # Kinetic: 0.6 along t0, 0.3 along t1
    
    print(f"  Anisotropic static friction: {mu_s_aniso}")
    print(f"  Anisotropic kinetic friction: {mu_k_aniso}")
    print(f"  (First component = friction along first tangent direction)")
    print(f"  (Second component = friction along second tangent direction)")
    
    # Assign to each collision
    for i in range(tangential_collisions.size()):
        tangential_collisions[i].mu_s_aniso = mu_s_aniso
        tangential_collisions[i].mu_k_aniso = mu_k_aniso
        print(f"  Assigned anisotropic coefficients to collision {i}")
    print()
    
    # -------------------------------------------------------------------------
    # Step 4: Demonstrate Effective Friction for Different Directions
    # -------------------------------------------------------------------------
    print("Step 4: Computing effective friction for different velocity directions...")
    print()
    
    # Test different velocity directions (angles from 0 to 2π)
    num_directions = 16
    angles = np.linspace(0, 2 * np.pi, num_directions, endpoint=False)
    
    print("Direction-dependent effective friction coefficients:")
    print("-" * 70)
    print(f"{'Angle (deg)':>12} {'Direction':>20} {'μ_s_eff':>12} {'μ_k_eff':>12}")
    print("-" * 70)
    
    mu_s_eff_values = []
    mu_k_eff_values = []
    
    for angle in angles:
        # Create unit direction vector in tangent plane
        # This represents the direction of tangential velocity
        tau_dir = np.array([np.cos(angle), np.sin(angle)])
        
        # Compute effective friction coefficients using elliptical model
        # μ_eff = sqrt((μ₀ t₀)² + (μ₁ t₁)²)
        mu_s_eff, mu_k_eff = ipctk.anisotropic_mu_eff_sqrt_mu0_t0_sq_plus_mu1_t1_sq(
            tau_dir, mu_s_aniso, mu_k_aniso
        )
        
        mu_s_eff_values.append(mu_s_eff)
        mu_k_eff_values.append(mu_k_eff)
        
        # Format direction vector for display
        dir_str = f"({tau_dir[0]:.3f}, {tau_dir[1]:.3f})"
        print(f"{np.degrees(angle):12.1f} {dir_str:>20} {mu_s_eff:12.4f} {mu_k_eff:12.4f}")
    
    print("-" * 70)
    print()
    
    # -------------------------------------------------------------------------
    # Step 5: Compare with Isotropic Case
    # -------------------------------------------------------------------------
    print("Step 5: Comparison with isotropic friction...")
    
    # For isotropic friction, effective mu is constant regardless of direction
    isotropic_mu_s = 0.5
    isotropic_mu_k = 0.3
    
    print(f"  Isotropic static friction: {isotropic_mu_s} (constant for all directions)")
    print(f"  Isotropic kinetic friction: {isotropic_mu_k} (constant for all directions)")
    print()
    print("  Anisotropic friction varies with direction:")
    print(f"    μ_s_eff range: [{min(mu_s_eff_values):.4f}, {max(mu_s_eff_values):.4f}]")
    print(f"    μ_k_eff range: [{min(mu_k_eff_values):.4f}, {max(mu_k_eff_values):.4f}]")
    print()
    
    # Show specific examples
    print("  Examples:")
    # Direction along first tangent (angle = 0)
    tau_dir_0 = np.array([1.0, 0.0])
    mu_s_eff_0, mu_k_eff_0 = ipctk.anisotropic_mu_eff_sqrt_mu0_t0_sq_plus_mu1_t1_sq(
        tau_dir_0, mu_s_aniso, mu_k_aniso
    )
    print(f"    Direction (1, 0): μ_s_eff = {mu_s_eff_0:.4f}, μ_k_eff = {mu_k_eff_0:.4f}")
    print(f"      (Should equal mu_s_aniso[0] = {mu_s_aniso[0]:.4f} and mu_k_aniso[0] = {mu_k_aniso[0]:.4f})")
    
    # Direction along second tangent (angle = π/2)
    tau_dir_90 = np.array([0.0, 1.0])
    mu_s_eff_90, mu_k_eff_90 = ipctk.anisotropic_mu_eff_sqrt_mu0_t0_sq_plus_mu1_t1_sq(
        tau_dir_90, mu_s_aniso, mu_k_aniso
    )
    print(f"    Direction (0, 1): μ_s_eff = {mu_s_eff_90:.4f}, μ_k_eff = {mu_k_eff_90:.4f}")
    print(f"      (Should equal mu_s_aniso[1] = {mu_s_aniso[1]:.4f} and mu_k_aniso[1] = {mu_k_aniso[1]:.4f})")
    
    # Diagonal direction (angle = π/4)
    tau_dir_45 = np.array([1.0, 1.0]) / np.sqrt(2.0)
    mu_s_eff_45, mu_k_eff_45 = ipctk.anisotropic_mu_eff_sqrt_mu0_t0_sq_plus_mu1_t1_sq(
        tau_dir_45, mu_s_aniso, mu_k_aniso
    )
    print(f"    Direction (1/√2, 1/√2): μ_s_eff = {mu_s_eff_45:.4f}, μ_k_eff = {mu_k_eff_45:.4f}")
    print(f"      (Should be sqrt(({mu_s_aniso[0]:.1f}·1/√2)² + ({mu_s_aniso[1]:.1f}·1/√2)²) = {np.sqrt((mu_s_aniso[0]/np.sqrt(2))**2 + (mu_s_aniso[1]/np.sqrt(2))**2):.4f})")
    print()
    
    # -------------------------------------------------------------------------
    # Step 6: Summary
    # -------------------------------------------------------------------------
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    print("Anisotropic friction allows friction to vary with the direction")
    print("of tangential velocity. The effective friction coefficient is")
    print("computed using the elliptical L2 projection model:")
    print()
    print("  μ_eff = sqrt((μ₀ t₀)² + (μ₁ t₁)²)")
    print()
    print("where t = τ / ||τ|| is the unit direction vector of tangential")
    print("velocity, and μ₀, μ₁ are the friction coefficients along the")
    print("two tangent basis directions.")
    print()
    print("This creates an elliptical friction model where friction is")
    print("strongest along the axes defined by mu_s_aniso and mu_k_aniso.")
    print("=" * 70)


if __name__ == "__main__":
    main()
