# Cookie Benchmark
This is an implementation of the [cookie bencmark](https://github.com/UM-Bridge/benchmarks/tree/main/benchmarks/cookies-problem-propagation) test problem.
This sets up a Docker container ```benmkent/cookiebenchmark:latest``` for use with the [UM-BRIDGE](https://um-bridge-benchmarks.readthedocs.io/en/docs/) uncertainty quantification interface.
This is easily downloaded and run using
``` docker run -p 4243:4242 -itd benmkent/cookiebenchmark:latest```
where the first port number is mapped through to ```4242``` which is exposed in the container.

The Docker container serves four models
- forward: the forward elliptic model with configurable model inputs
- benchmark: the forward elliptic model with fixed model inputs
- cookietime: the forward parabolic model with configurable inputs
- cookietimebenchmark: the forward parabolic model with fixed model inputs

# Implementation details
## ellipticpde.py
A Python implementation of the elliptic and parabolic ``cookie'' PDE problem.
Finite element approximation on a quad grid is assembled using FEniCS.
The linear systems are solved using PETSc via the petsc4py interface.

## umbridge-server.py
This defines the UM-BRIDGE interfaces for the test problem.

Model |Dictionary Key | Default value | User control | Description
------|---------------|---------------|--------------|-------------
Both | N | 400 | 100 * config['Fidelity'] or config['N'] integer | The number of cells in each dimension (i.e. a mesh of N^2 elements)
Both | BasisDegree | 1 | Integer | The degree of the piecewise polynomial FE approximation
Both | quad_degree | 8 | Integer | The quadrature degree used to evaluate the forcing and coefficient functions
Both | coeffs | None | Float Vector | coeff[0] defines the background diffusion field (equal to 1.0 bu default)
Elliptic | pc  | "none" | "ILU" or "JACOBI" | Preconditioning for the GM-RES solver
Elliptic | tol | "LU" | Float | Relative tolerance for the GM-RES solver. "LU" uses exact solution via full LU preconditioning.
Parabolic | letol  | 1e-4 | Float | Local timestepping error tolerance for a simple implementation of TR-AB2 timestepping.
Parabolic | T | 10.0 | Float | Final time for timestepping approximation. The QoI is evaluated and returned for time T.

## run_forward_benchmark_in_matlab_fenics.m
Runs the FEniCS implementation of the benchmark problem in MATLAB via the UM-BRIDGE interface and a Docker running ```benmkent/cookiebenchmark:latest```.
Requires the [Sparse Grids MATLAB Kit](https://sites.google.com/view/sparse-grids-kit).

## test_output.py
Python script to generate QoI approximations for testing.
