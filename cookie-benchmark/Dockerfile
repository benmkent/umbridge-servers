FROM ubuntu:latest

# Install fenics and python packages
RUN apt update
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:fenics-packages/fenics
RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt install -y bash fenics python3-pip
RUN pip3 install --break-system-packages umbridge

# Expose port for UM-BRIDGE server
EXPOSE 4242

RUN apt update

# Set up internal folders and files
WORKDIR /
COPY umbridge_server.py /
COPY cookiepde.py /

# Start Python UM-BRDIGE server
CMD python3 umbridge_server.py
