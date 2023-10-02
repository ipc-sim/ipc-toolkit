---
title: 'IPC Toolkit: A Modular Library for Intersection-free Contact Handling'
tags:
  - C++
  - Python
  - physical simulation
  - frictional contact
  - collision detection
authors:
  - name: Zachary Ferguson
    orcid: 0000-0003-2466-3768
    affiliation: "1, 2"
affiliations:
 - name: New York University
   index: 1
 - name: Massachusetts Institute of Technology
   index: 2
date: October 1, 2023
bibliography: ./source/refs.bib
---

# Summary

Physical simulation has become a critical component of many scientific and engineering applications. In particular, the simulation of frictional contact is a key component of many applications, including computer graphics, robotics, biomechanics, and computational physics. However, the simulation of frictional contact is a challenging problem, and many existing methods allow for interpenetration of objects, which can lead to a variety of issues, including non-physical behavior, numerical instability, and even simulation failure.

Incremental Potential Contact (IPC) [@Li2020IPC] robustly handles frictional contact without allowing for interpenetration.

# Statement of need

We integrate techniques from [@Li2020IPC; @Li2021CIPC; @Li2023Convergent; @Ferguson2023HighOrderIPC]. We also utilize the provably conservative continuous collision detection (CCD) methods of @Wang2021TightInclusion and @Belgrod2023Time.

The IPC Toolkit has been used in a number of research projects [@Ferguson2021RigidIPC; @Huang2022Differentiable; @Ferguson2023HighOrderIPC; @Ferguson2023InTimestepRemeshing; @Vidulis2023Computational; @Tang2023Chainmail], and is currently utilized in the PolyFEM [@polyfem] and Rigid IPC [@Ferguson2021RigidIPC] libraries.

# Acknowledgements

This work was also partially supported by the NSF CAREER award under Grant No. 1652515, the NSF grants OAC-1835712, OIA-1937043, CHS-1908767, CHS-1901091, NSERC DGECR-2021-00461 and RGPIN 2021-03707, a Sloan Fellowship, a gift from Adobe Research and a gift from Advanced Micro Devices, Inc.

# References