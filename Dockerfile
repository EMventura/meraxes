# Install required tools
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /meraxes

# Copy all source code into container
COPY . .

# Create and move into build dir
RUN mkdir build && cd build && cmake .. && make -j$(nproc)

# Set default command (meraxes is the output binary file name)
CMD ["./build/meraxes"]
