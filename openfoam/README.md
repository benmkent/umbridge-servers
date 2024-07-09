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
The compressed size of the container is 789.48 MB.

The Docker container a single models
- forward

## Properties
The parametric input is two dimensional.
The parameters are
- jet flow rate, values in [2.34, 23.4]
- inflow rate, values in [3.46, 34.6].
The output is the reattachement point of the flow.

Mapping | Dimensions	| Description
--------|-------------|------------
input |	[2] |	These values modify the flow properties (software does not check that input values are valid).
output |	[1] |	The reattachement point of the flow after the hump.

The model supports the UM-BRIDGE evaluate feature.

Feature	| Supported
--------|---------
Evaluate|	True
Gradient|	False
ApplyJacobian|	False
ApplyHessian|	False

### Config
Dictionary Key | Default value | User control | Description
---------------|---------------|--------------|-------------
Fidelity       | n/a           | 1,2,3,4      | Selects one of four meshes, with coarsest fidelity 1 and finest fidelity 4.
res_tol        |  1e-10        | Float        | Tol on the residual in the iterative solvers in the simpleFoam solver.

## Implementation details
### umbridge-server.py
This defines the UM-BRIDGE interfaces for the test problem.
This copies the case file for the selected fidelity into a temporary folder.
The simpleFoam solver is applied to the temporary folder.
5000 iterations are performed.
The wallShearStress is extracted and posprocessed using Python.
It is saved to the directory ```./outputdata```.
A cubic spline is used to interpolate the discrete wallShearStress data along the bottomWall of the domain.
The final root (change in sign of wallShearStress) indicates the reattachment point.
This is the QoI and is returned.

### NASA_hump_data_*
These files contain the blockMesh definitions for the four fidelities.
The meshes are constructed and an initial condition prescribed by the data in NASA_hump_data_baseline.
This occurs during the building of the Docker image (see Dockerfile).

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
