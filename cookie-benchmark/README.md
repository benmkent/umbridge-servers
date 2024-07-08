# Cookie Benchmark
## Overview
This is an implementation of the [cookie bencmark](https://github.com/UM-Bridge/benchmarks/tree/main/benchmarks/cookies-problem-propagation) test problem.
This sets up a Docker container ```benmkent/cookiebenchmark:latest``` for use with the [UM-BRIDGE](https://um-bridge-benchmarks.readthedocs.io/en/docs/) uncertainty quantification interface.
This is easily downloaded and run using
``` docker run -p 4242:4242 -itd benmkent/cookiebenchmark:latest```
where the first port number is mapped through to ```4242``` which is exposed in the container.

## Run 
```
docker run -p 4242:4242 -it benmkent/cookiebenchmark:latest
```
The compressed size of the container is 754.5 MB.

The Docker container serves four models
- forward: the forward elliptic model with configurable model inputs
- benchmark: the forward elliptic model with fixed model inputs
- cookietime: the forward parabolic model with configurable inputs
- cookietimebenchmark: the forward parabolic model with fixed model inputs

## Properties
All four models take the same parametric input and return the same QoI (the integral of the solution over the spatial subdomain [0.4,0.6]^2.

Mapping | Dimensions	| Description
--------|-------------|------------
input |	[8] |	These values modify the conductivity coefficient in the 8 cookies, each of them must be greater than -1 (software does not check that input values are valid)
output |	[1] |	The integral of the solution over the central subdomain

All four models suppor the UM-BRDIGE evaluate feature.

Feature	| Supported
--------|---------
Evaluate|	True
Gradient|	False
ApplyJacobian|	False
ApplyHessian|	False

### Config
Model |Dictionary Key | Default value | User control | Description
------|---------------|---------------|--------------|-------------
Both | N | 400 | 100 * config['Fidelity'] or config['N'] integer | The number of cells in each dimension (i.e. a mesh of N^2 elements)
Both | BasisDegree | 1 | Integer | The degree of the piecewise polynomial FE approximation
Both | quad_degree | 8 | Integer | The quadrature degree used to evaluate integrals in the matrix assembly.
Both | coeffs | None | Float Vector | coeff[0] defines the background diffusion field (equal to 1.0 bu default)
Elliptic | pc  | "none" | "ILU" or "JACOBI" | Preconditioning for the GM-RES solver
Elliptic | tol | "LU" | Float | Relative tolerance for the GM-RES solver. "LU" uses exact solution via full LU preconditioning.
Parabolic | letol  | 1e-4 | Float | Local timestepping error tolerance for a simple implementation of TR-AB2 timestepping.
Parabolic | T | 10.0 | Float | Final time for timestepping approximation. The QoI is evaluated and returned for time T.

## Implementation details
### ellipticpde.py
A Python implementation of the elliptic and parabolic ``cookie'' PDE problem.
Finite element approximation on a quad grid is assembled using FEniCS.
The linear systems are solved using PETSc via the petsc4py interface.

### umbridge-server.py
This defines the UM-BRIDGE interfaces for the test problem.

### run_forward_benchmark_in_matlab_fenics.m
Runs the FEniCS implementation of the benchmark problem in MATLAB via the UM-BRIDGE interface and a Docker running ```benmkent/cookiebenchmark:latest```.
Requires the [Sparse Grids MATLAB Kit](https://sites.google.com/view/sparse-grids-kit).

### test_output.py
Python script to generate QoI approximations for testing. This evaluates the quantity of interest for three different parameters at four different fidelities.
Results can be piped from the console for plotting
```
python3 test_output.py http://0.0.0.0:4242 >> results.txt
```