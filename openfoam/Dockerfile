FROM ubuntu:mantic

RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt install -y bash python3-pip git curl wget
# Add the repository
RUN curl https://dl.openfoam.com/add-debian-repo.sh | bash
RUN apt-get install -y openfoam2312-default

RUN pip3 install --break-system-packages umbridge

EXPOSE 4242

RUN apt update

WORKDIR /

# Copy files and subdirectories (i.e. the openfoam setup)
COPY . .
RUN wget -P ./2DWMH/pitzDailyModified https://turbmodels.larc.nasa.gov/Nasahump_grids/hump2newtop_noplenumZ409x109.nmf   

RUN wget -P ./2DWMH/pitzDailyModified https://turbmodels.larc.nasa.gov/Nasahump_grids/hump2newtop_noplenumZ409x109.p3dfmt.gz 
RUN gzip -dfk ./2DWMH/pitzDailyModified/hump2newtop_noplenumZ409x109.p3dfmt.gz

RUN python3 ./p3d2gmsh.py -o ./2DWMH/pitzDailyModified/hump2newtop_noplenumZ409x109.msh -m hump2newtop_noplenumZ409x109.nmf hump2newtop_noplenumZ409x109.p3dfmt
RUN openfoam2312 gmshToFoam -case ./2DWMH/pitzDailyModified/ hump2newtop_noplenumZ409x109.msh

# Docker file should pull meshes from e.g. https://turbmodels.larc.nasa.gov/Nasahump_grids/hump2newtop_noplenumZ1633x433.p3dfmt.gz and convert

CMD python3 umbridge-server.py