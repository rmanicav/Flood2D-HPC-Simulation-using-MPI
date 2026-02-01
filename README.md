# Flood2D HPC Simulation using MPI

## Overview
Parallel C++ simulation of 2D flood propagation using MPI.

Flood2D is a research-oriented scientific computing project implementing a
two-dimensional flood simulation using MPI-based parallelization in C++.

## Key Features
- 2D grid-based numerical simulation
- Domain decomposition across MPI processes
- Parallel execution using MPI (message passing)
- Deterministic and reproducible simulation results
- Designed for execution on HPC clusters and multi-core systems

## Technical Approach
- C++ implementation with explicit memory management
- MPI for inter-process communication and synchronization
- Structured grid decomposition across ranks
- Boundary data exchange between neighboring subdomains
- Scalable execution model suitable for distributed environments

## Technologies
- **Language:** C++
- **Parallel Computing:** MPI
- **Execution Environment:** Linux, HPC clusters
- **Build & Run:** Standard MPI toolchain (`mpic++`, `mpirun`)

## Project Structure
src/        -> source code  
results/    -> outputs  
assets/     -> figures / performance plots  


## Example Execution
mpic++ flood2d.cpp -o flood2d
mpirun -np 4 ./flood2d


## Research Software Engineering Context

This project was developed to demonstrate parallel scientific software
design, HPC execution models, and reproducible research code.

It serves as a reference example for Research Software Engineer and
scientific computing roles where scalable computation and collaboration
with domain scientists are essential.

## Skills Demonstrated
- C++ scientific programming
- MPI-based parallel computing
- Domain decomposition
- HPC execution workflows
- Research software engineering

## Authorship & Contributions
This project was developed in a collaborative research environment.
The overall simulation framework and architectural design were created
as part of a broader research effort, while I contributed to the
implementation, MPI-based parallelization, integration, execution
workflow, and research-oriented adaptation of the codebase.

This repository is shared to demonstrate experience with scientific
computing, MPI-based parallel programming, and research software
engineering practices.



Dr. Rajesh Manicavasagam
Research Software Engineer | Scientific Computing | HPC & MPI
GitHub: https://github.com/rmanicav
