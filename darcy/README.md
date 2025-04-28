# Darcy Flow
## Overview
An UM-BRIDGE server can be run using Docker with the command
``` docker run -p 4242:4242 -itd benmkent/darcy:latest```
where the first port number is mapped through to ```4242``` which is exposed in the container.

## Authors
- [Benjamin Kent](mailto:kent@imati.cnr.it)

## Run 
```
docker run -p 4242:4242 -it benmkent/darcy:latest
```
The compressed size of the container is not recorded yet.

The Docker container serves two models.
- `darcy`
- `advdiff`

## Properties

Mapping | Dimensions	| Description
--------|-------------|------------
input |	[1,50,50,50] |	
output |	[50] |	

The models only support the UM-BRDIGE evaluate feature.

Feature	| Supported
--------|---------
Evaluate|	True
Gradient|	False
ApplyJacobian|	False
ApplyHessian|	False

### Config
The basic configuration parameters are:

Model |Dictionary Key | Default value | User control | Description
------|---------------|---------------|--------------|-------------


## Mount directories
Mount directory | Purpose
---             |---
None            | 

## Description

## Implementation details
