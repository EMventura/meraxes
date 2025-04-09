FROM ubuntu:22.04

# Install required tools and MPI
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    mpich \
    libmpich-dev \
    pkg-config \
    libgsl-dev \
    libfftw3-dev \
    libfftw3-mpi-dev \
    libhdf5-mpi-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /meraxes

# Copy all source code into container
COPY . .

# Configure and build
RUN mkdir build && cd build && cmake .. && make -j$(nproc)

# Set default command
CMD ["./build/meraxes"]

