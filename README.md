# Cookie Benchmark
This is an implementation of the [cookie bencmark](https://github.com/UM-Bridge/benchmarks/tree/main/benchmarks/cookies-problem-propagation) test problem.
This sets up a Docker container ```benmkent/cookiebenchmark:latest``` for use with the [UM-BRIDGE](https://um-bridge-benchmarks.readthedocs.io/en/docs/) uncertainty quantification interface.

# Implementation details
## ellipticpde.py
A Python implementation of the elliptic PDE problem.
Piecewise linear approximation finite element approximation is assembled using FEniCS.
The linear systems are solved using PETSc via the petsc4py interface.

## umbridge-server.py
This defines the UM-BRIDGE interfaces for the test problem.
