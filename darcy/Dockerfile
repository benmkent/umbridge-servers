FROM ubuntu:latest

# Set environment variables to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies for Julia, Python, and other necessary packages
RUN apt-get update && \
    apt-get install -y \
    wget \
    ca-certificates \
    build-essential \
    git \
    liblapack-dev \
    libopenblas-dev \
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    libpython3-dev \
    curl \
    && apt-get clean

# Install Julia
RUN curl -fsSL https://install.julialang.org | sh -s -- -y

# Add Julia to PATH
ENV PATH="/root/.juliaup/bin:${PATH}"

# Set the environment variable for OpenSSL version (ignore or specify version)
ENV JULIA_CONDAPKG_OPENSSL_VERSION="ignore"


# Set up internal folders and files
WORKDIR /
COPY . .

RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile();'

# Start Julia server
CMD julia -e 'using Pkg; Pkg.activate("."); include("DarcyServer.jl");'