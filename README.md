[![Build and Push Docker Image (OpenFOAM)](https://github.com/benmkent/umbridge-servers/actions/workflows/docker-image-openfoam.yml/badge.svg)](https://github.com/benmkent/umbridge-servers/actions/workflows/docker-image-openfoam.yml)
[![Build and Push Docker Image (Cookies)](https://github.com/benmkent/umbridge-servers/actions/workflows/docker-image-cookie.yml/badge.svg)](https://github.com/benmkent/umbridge-servers/actions/workflows/docker-image-cookie.yml)
[![Build and Push Docker Image (Double Glazing)](https://github.com/benmkent/umbridge-servers/actions/workflows/docker-image-doubleglazing.yml/badge.svg)](https://github.com/benmkent/umbridge-servers/actions/workflows/docker-image-doubleglazing.yml)
[![Build and Push Docker Image (Darcy)](https://github.com/benmkent/umbridge-servers/actions/workflows/darcy.yml/badge.svg)](https://github.com/benmkent/umbridge-servers/actions/workflows/darcy.yml)

# UM-BRIDGE Servers

This repository contains the code and Dockerfiles for containers to be used with [UM-BRIDGE](https://github.com/UM-Bridge/umbridge).

## Setup
In the simplest instance, these servers can be used as described in the [UM-BRIDGE quickstart guide](https://um-bridge-benchmarks.readthedocs.io/en/docs/quickstart.html).

To set up to use a Kubernetes cluster the following steps are recommended.
The simplest set up follows the recommenedations [UM-BRIDGE Kubernetes](https://um-bridge-benchmarks.readthedocs.io/en/docs/umbridge/kubernetes.html)

1. Install and set up a Kubernetes backend e.g. [micro-k8s](https://microk8s.io/docs/getting-started)
2. Navigate to the model folder e.g. `cd ./doubleglazing'
3. (Ignore for `cookies`) Create the persistent volume and persistent volume claim. This allows the Kubernetest pods to read and write data in the `/output_data` subfolders.
   - `microk8s kubectl create -f pv.yaml`
   - `microk8s kubectl create -f pvc.yaml`
4. Create the pods for the model.
   - `microk8s kubectl create -f model.yaml`
   - The file `model.yaml' allows the number of pods, the resources etc to be tweaked.
5. To delete use `microk8s kubectl delete -f model.yaml`.

## Setup of `svc-ingress.yaml`
To run multiple models through the same microk8s set up and IP address and ports we use the following set up.

There are two components for each model in the `svc-ingress.yaml` file.
Firstly services are set up for a model, in this case `openfoam`.
This selects `app: openfoam`.
```
# Service for OpenFOAM application
apiVersion: v1
kind: Service
metadata:
  name: svc-openfoam
  namespace: default
  annotations:
    haproxy.org/load-balance: "leastconn"  # Use least connections for load balancing
spec:
  ports:
    - port: 80
      protocol: TCP
      targetPort: 4242
  selector:
    app: openfoam
```
Ingress specifies how this service is reached.
The ingress has a `rule` identifying paths `/openfoam` that connect to the Kubernetes access port.
The `http` message is passed on, with `/openfoam` stripped from the message so that it can be properly parsed by UM-BRIDGE.
```
# Ingress for OpenFOAM application
apiVersion: networking.k8s.io/v1
kind: Ingress
metadata:
  name: openfoam-ingress
  annotations:
    haproxy.org/path-rewrite: /openfoam/(.*) /\1
spec:
  ingressClassName: haproxy
  rules:
    - http:
        paths:
          - path: /openfoam
            pathType: Prefix
            backend:
              service:
                name: svc-openfoam
                port:
                  number: 80
```
