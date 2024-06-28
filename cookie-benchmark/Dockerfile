FROM ubuntu:mantic

RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt install -y bash fenics python3-pip
RUN pip3 install --break-system-packages umbridge

EXPOSE 4242

RUN apt update

WORKDIR /

COPY umbridge-server.py /
COPY ellipticpde.py /

CMD python3 umbridge-server.py