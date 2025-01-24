# Openfoam NASA 2DWMH Flow Seperation
## Overview
This is an implementation of the [NASA 2D wall mounted hump](https://github.com/UM-Bridge/benchmarks/tree/main/benchmarks/cookies-problem-propagation](https://turbmodels.larc.nasa.gov/nasahump_val.html) test problem.
This sets up a Docker container ```benmkent/openfoam2dwmh:latest``` for use with the [UM-BRIDGE](https://um-bridge-benchmarks.readthedocs.io/en/docs/) uncertainty quantification interface.
This is easily downloaded and run using
``` docker run -p 4242:4242 -itd benmkent/openfoam2dwmh:latest```
where the first port number is mapped through to ```4242``` which is exposed in the container.

## Run 
```
docker run -p 4242:4242 -it benmkent/openfoam2dwmh:latest
```
The compressed size of the container is 787.37 MB.

The Docker container a single model
- forward2dwmh

## Properties
The parametric input is two dimensional.
The parameters are
- jet flow rate,
- inflow rate.
A number of quantities of interest can be evaluated on the solution.

Mapping | Dimensions	| Description
--------|-------------|------------
input |	[2] |	These values modify the flow properties (software does not check that input values are valid).
output |	[varied] |	Various quantites of interest as detailed below.

The model supports the UM-BRIDGE evaluate feature.

Feature	| Supported
--------|---------
Evaluate|	True
Gradient|	False
ApplyJacobian|	False
ApplyHessian|	False

### Config
Dictionary Key | Default value | User control          | Description
---------------|---------------|-----------------------|-------------
Fidelity       |               | 20                    | Fine mesh (January 25)
               |               | 11                    | Intermediate mesh (January 25)
               |               | 12                    | Coarse mesh (January 25)
               |               | 4                     | Baseline mesh
qoi            |               | "reattachmentpoint"   | Flow reattachment point
               |               | "exectime"            | Solver execution time
               |               | "cf"                  | Friction coefficient on bottom boundary
               |               | "cp"                  | Pressure coefficient on bottom boundary
               |               | "yplus"               | yPlus on bottom boundary
               |               | "p"                   | Pressure on bottom boundary
               |               | "forces"              | Forces on hump (total, pressure, viscous)
               |               | "residuals"           | Solver residuals at each simplefoam iteration
               |               | "yslice"              | Ux at a vertical slice in the domain (default x=0.5, otherwise specify "yslice_x" in config.
abs_tol        |  1e-10         | Float                | Absolute tolerance for iterative LA solvers in each simplefoam iteration.
rel_tol        |  1e-3         | Float                 | Relative tolerance for iterative LA solvers in each simplefoam iteration.
res_tol_u      |  1e-10         | Float                | Tolerance for velocity residuals in simplefoam.
res_tol_p      |  1e-10         | Float                | Tolerance for pressure residual in simplefoam.
res_tol_nut    |  1e-10         | Float                | Tolerance for nuTilda residual in simplefoam.

Notes:
- Unintuative mesh numbering is a consequence of legacy meshes.
- Computation of forces QoI is not availible on the baseline mesh.

## Implementation details
### umbridge-server.py
This defines the UM-BRIDGE interfaces for the test problem.
This copies the case file for the selected fidelity into a temporary folder.
The simpleFoam solver is applied to the temporary folder.
Either 5000 iterations are performed or the case terminates early if the residaul tolerances are satisfied.
OpenFOAM postprocessing is applied.
It is saved to the directory ```./outputdata```.
Python based post processing extracts the QoI.

### NASA_hump_data_*
These files contain the case definitions for each of the fidelities.

### Kubernetes files
pv.yaml and pvc.yaml are used to bind the directory ```./outputdata``` to a folder on the host machine.
model.yaml creates the Kubernetes pods for parallel processing.
To start the cluster
```
kubectl apply -f pv.yaml
kubectl apply -f pvc.yaml
kubectl apply -f model.yaml
```
This is then deleted using ```kubectl delete -f model.yaml```.
