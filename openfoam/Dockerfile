FROM ubuntu:mantic

RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt install -y bash python3-pip git curl wget gmsh
# Add the repository
RUN curl https://dl.openfoam.com/add-debian-repo.sh | bash
RUN apt-get install -y openfoam2406-dev

RUN pip3 install --break-system-packages umbridge numpy fluidfoam

EXPOSE 4242

RUN apt update

WORKDIR /


# Copy files and subdirectories (i.e. the openfoam setup)
COPY . .

# Make sure we build the jetNasaHump tool
RUN openfoam2406 wclean ./jetNasaHump
RUN openfoam2406 wmake ./jetNasaHump

# Map fields from common initial guess and ensure meshes are OK
RUN openfoam2406 blockMesh -case ./NASA_hump_data_baseline

RUN openfoam2406 blockMesh -case ./NASA_hump_data_coarse1
RUN openfoam2406 mapFieldsPar -consistent -case ./NASA_hump_data_coarse1 ./NASA_hump_data_baseline
RUN openfoam2406 blockMesh -case ./NASA_hump_data_coarse2
RUN openfoam2406 mapFieldsPar -consistent -case ./NASA_hump_data_coarse2 ./NASA_hump_data_baseline
RUN openfoam2406 blockMesh -case ./NASA_hump_data_coarse3
RUN openfoam2406 mapFieldsPar -consistent -case ./NASA_hump_data_coarse3 ./NASA_hump_data_baseline

# Once openfoam folders are configured we can start the server
CMD python3 umbridge-server.py