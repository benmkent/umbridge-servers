# Cookie Benchmark
This is an implementation of the [cookie bencmark](https://github.com/UM-Bridge/benchmarks/tree/main/benchmarks/cookies-problem-propagation) test problem.
This sets up a Docker container ```benmkent/cookiebenchmark:latest``` for use with the [UM-BRIDGE](https://um-bridge-benchmarks.readthedocs.io/en/docs/) uncertainty quantification interface.
This is easily downloaded and run using
``` docker run -p 4243:4242 -itd benmkent/cookiebenchmark:latest```
where the first port number is mapped through to ```4242``` which is exposed in the container.

# Implementation details
## ellipticpde.py
A Python implementation of the elliptic PDE problem.
Piecewise linear approximation finite element approximation is assembled using FEniCS.
The linear systems are solved using PETSc via the petsc4py interface.

## umbridge-server.py
This defines the UM-BRIDGE interfaces for the test problem.

## run_forward_benchmark_in_matlab_fenics.m
Runs the FEniCS implementation of the benchmark problem in MATLAB via the UM-BRIDGE interface and a Docker running ```benmkent/cookiebenchmark:latest```.
Requires the [Sparse Grids MATLAB Kit](https://sites.google.com/view/sparse-grids-kit).
