# Flood2D — Parallel 2D Flood Simulation using C++ and MPI

## Overview
Flood2D is a research-oriented **scientific computing** project implementing a
two-dimensional flood simulation using **MPI-based parallelization** in C++.
The project focuses on **distributed numerical computation**, scalability, and
reproducibility in **HPC environments**.

The primary goal is to demonstrate **research software engineering practices**
for parallel simulations rather than domain-specific flood modeling novelty.

---

## Key Features
- 2D grid-based numerical simulation
- Domain decomposition across MPI processes
- Parallel execution using MPI (message passing)
- Deterministic and reproducible simulation results
- Designed for execution on HPC clusters and multi-core systems

---

## Technical Approach
- C++ implementation with explicit memory management
- MPI for inter-process communication and synchronization
- Structured grid decomposition across ranks
- Boundary data exchange between neighboring subdomains
- Scalable execution model suitable for distributed environments

---

## Technologies
- **Language:** C++
- **Parallel Computing:** MPI
- **Execution Environment:** Linux, HPC clusters
- **Build & Run:** Standard MPI toolchain (`mpic++`, `mpirun`)

---

## Example Execution
  bash
  mpirun -np 4 ./flood2d 

## Research Software Engineering Context

This project was developed to demonstrate parallel scientific software
design, HPC execution models, and reproducible research code.

It serves as a reference example for Research Software Engineer and
scientific computing roles where scalable computation and collaboration
with domain scientists are essential.

Author

Dr. Rajesh Manicavasagam
Research Software Engineer | Scientific Computing | HPC & MPI
GitHub: https://github.com/rmanicav
