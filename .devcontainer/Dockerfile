FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV CMAKE_VERSION=3.31.0

# Basic dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    wget \
    curl \
    ninja-build \
    ccache \
    libeigen3-dev \
    libtbb-dev \
    libspdlog-dev \
    python3 \
    python3-pip \
    python3-dev \
    libgmp-dev \
    libssl-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libxml2-dev \
    libjsoncpp-dev \
    libz3-dev \
    sudo \
    && rm -rf /var/lib/apt/lists/*

# Create a new user with sudo privileges
RUN useradd -m devuser && echo "devuser:password" | chpasswd && usermod -aG sudo devuser

# Configure passwordless sudo for convenience
RUN echo "devuser ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# Set up Python tools
RUN pip3 install --upgrade pip setuptools wheel pre-commit

ENV CCACHE_DIR=/root/.ccache
ENV CCACHE_MAXSIZE=1G
RUN mkdir -p $CCACHE_DIR

ENV CXX_STANDARD=17

# Install CMake
RUN mkdir -p /root/src && \
    cd /root/src && \
    wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}.tar.gz && \
    tar -xzvf cmake-${CMAKE_VERSION}.tar.gz && \
    cd cmake-${CMAKE_VERSION} && \
    ./bootstrap --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd /root && \
    rm -rf /root/src/cmake-${CMAKE_VERSION} /root/src/cmake-${CMAKE_VERSION}.tar.gz

RUN cmake --version

# Install LLVM/Clang and Clang-Format version 18
RUN wget -qO- https://apt.llvm.org/llvm.sh | sudo bash -s -- 18 && \
    apt-get install -y clang-format-18

# Create an alias for clang-format to clang-format-18
RUN ln -s /usr/bin/clang-format-18 /usr/bin/clang-format

RUN clang-18 --version && clang-format --version

# Cleanup to reduce image size
RUN apt-get remove -y build-essential git wget curl && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /root/src

# Set the default user and working directory
USER devuser
WORKDIR /home/devuser

CMD ["bash"]
