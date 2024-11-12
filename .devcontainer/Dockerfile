FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV CMAKE_VERSION=3.31.0
ENV LLVM_VERSION=18.0.0

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

RUN pip3 install --upgrade pip setuptools wheel pre-commit

ENV CCACHE_DIR=/root/.ccache
ENV CCACHE_MAXSIZE=1G
RUN mkdir -p $CCACHE_DIR

ENV CXX_STANDARD=17

WORKDIR /root

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

# Install LLVM
RUN mkdir -p /root/src && \
    cd /root/src && \
    wget https://github.com/llvm/llvm-project/releases/download/llvmorg-${LLVM_VERSION}/llvm-project-${LLVM_VERSION}.src.tar.xz && \
    tar -xf llvm-project-${LLVM_VERSION}.src.tar.xz && \
    mkdir -p llvm-project-${LLVM_VERSION}.src/build && \
    cd llvm-project-${LLVM_VERSION}.src/build && \
    cmake -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DLLVM_ENABLE_PROJECTS=clang -DCMAKE_INSTALL_PREFIX=/usr/local ../llvm && \
    ninja && \
    ninja install && \
    cd /root && \
    rm -rf /root/src/llvm-project-${LLVM_VERSION}.src /root/src/llvm-project-${LLVM_VERSION}.src.tar.xz

RUN clang++ --version
RUN clang-format --version

# Cleanup to reduce image size
RUN apt-get remove -y build-essential git wget curl && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /root/src

CMD ["bash"]
