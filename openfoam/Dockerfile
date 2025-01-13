FROM ubuntu:latest

RUN apt update

RUN DEBIAN_FRONTEND=noninteractive apt install -y bash python3-pip git openssl curl wget gmsh gnupg ca-certificates

# Add the repository
RUN apt update
# RUN curl https://dl.openfoam.com/add-debian-repo.sh
# RUN curl https://dl.openfoam.com/add-debian-repo.sh | bash
# RUN apt-get install -y openfoam2406-dev

WORKDIR /
# Copy files and subdirectories (i.e. the openfoam setup)
COPY . .

RUN bash addopenfoamrepo.sh
RUN apt update
RUN apt-get install -y openfoam2406-dev

# Set up python environment
RUN pip3 install --break-system-packages umbridge numpy fluidfoam
# Expose port for umbridge
EXPOSE 4242

RUN apt update

# Make sure we build the jetNasaHump tool
RUN openfoam2406 wclean ./jetNasaHump
RUN openfoam2406 wmake ./jetNasaHump

# Map fields from common initial guess and ensure meshes are OK
RUN openfoam2406 blockMesh -case ./NASA_hump_data_baseline

# Build Jan25 test meshes 
# This is currently zero - an initial guess, perhaps from a lofi solve could be used.
# These meshes seperate out the `hump' as a separate patch.
RUN openfoam2406 blockMesh -case ./NASA_hump_data_jan25_fine
#RUN openfoam2406 mapFieldsPar -consistent -case ./NASA_hump_data_jan25_fine ./NASA_hump_data_baseline

RUN openfoam2406 blockMesh -case ./NASA_hump_data_jan25_coarse
RUN openfoam2406 mapFieldsPar -consistent -case ./NASA_hump_data_jan25_coarse ./NASA_hump_data_jan25_fine

RUN openfoam2406 blockMesh -case ./NASA_hump_data_jan25_coarse_2
RUN openfoam2406 mapFieldsPar -consistent -case ./NASA_hump_data_jan25_coarse_2 ./NASA_hump_data_jan25_fine

# Once openfoam folders are configured we can start the server
CMD python3 umbridge_server.py