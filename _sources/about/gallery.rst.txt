Gallery
=======

The IPC Toolkit is used by a variety of projects across academia and industry for collision handling in physical simulations. If your project uses the IPC Toolkit, please `let us know <https://github.com/ipc-sim/ipc-toolkit/discussions>`_ and we'll add it here!

.. raw:: html
   :file: ../_static/gallery.html

Simulation Frameworks
---------------------

.. list-table::
   :widths: 25 75
   :header-rows: 0

   * - `PolyFEM <https://polyfem.github.io/>`_
     - A polyvalent C++ and Python finite element library supporting a wide range of PDEs including elasticity, Stokes, and Navier-Stokes. PolyFEM uses the IPC Toolkit for contact handling in its simulations.

   * - `SfePy <https://github.com/sfepy/sfepy>`_
     - Simple Finite Elements in Python — a software for solving systems of coupled partial differential equations by the finite element method in 1D, 2D, and 3D. It can be used both as a black-box PDE solver and as a Python package for building custom applications.

   * - `FEDOO <https://3mah.github.io/fedoo-docs>`_
     - A free open-source Python finite element library mainly dedicated to mechanical problems, including support for finite strain, composites, elasto-plastic laws, homogenization, and 2D/3D contact with self-contact.

   * - `Rigid IPC <https://github.com/ipc-sim/rigid-ipc>`_
     - Robust, intersection-free simulations of rigid bodies. This is the reference implementation of the SIGGRAPH 2021 paper *Intersection-free Rigid Body Dynamics*, and uses the IPC Toolkit for common IPC functions.

Papers Using IPC Toolkit
------------------------

The following is a non-exhaustive list of publications that use or build upon the IPC Toolkit:

- **Intersection-free Rigid Body Dynamics** :cite:`Ferguson2021RigidIPC`
- **High-Order Incremental Potential Contact for Elastodynamic Simulation on Curved Meshes** :cite:`Ferguson2023HighOrderIPC`
- **In-Timestep Remeshing for Contacting Elastodynamics** :cite:`Ferguson2023ITR`
- **Efficient Incremental Potential Contact for Actuated Face Simulation** :cite:`Li2023ActuatedFaceSimulation`
- **Computational Exploration of Multistable Elastic Knots** :cite:`Vidulis2023MultistableElasticKnots`
- **Differentiable solver for time-dependent deformation problems with contact** :cite:`Huang2024Differentiable`
- **A systematic comparison between FEBio and PolyFEM for biomechanical systems** :cite:`Martin2024SystematicComparison`
- **Geometric Contact Potential** :cite:`Huang2025GCP`
- **Adopting Research in Production at Netflix Animation Studios** :cite:`Andrus2025Adopting`
- **ANIME-Rod: Adjustable Nonlinear Isotropic Materials for Elastic Rods** :cite:`Chen2025ANIMERod`
- **Intersection-Free Garment Retargeting** :cite:`Huang2025Garment`

Add Your Project
----------------

If you are using the IPC Toolkit in your research or software, we would love to feature your work here. Please open a `discussion post <https://github.com/ipc-sim/ipc-toolkit/discussions>`_ or submit a `pull request <https://github.com/ipc-sim/ipc-toolkit/pulls>`_ adding your project to this page.