# CUDA-enabled development container for the IPC Toolkit.
#
# This mirrors ../Dockerfile (the default Ubuntu dev container) but is based on
# NVIDIA's CUDA "devel" image so nvcc is available. Like the default container
# it does NOT copy the source: the workspace is bind-mounted at runtime by the
# Dev Containers tooling. Compiling CUDA does not require a GPU; running it
# does. On a Linux host with the NVIDIA Container Toolkit you can expose the GPU
# by adding "--gpus=all" to runArgs in devcontainer.json.

ARG CUDA_IMAGE=nvidia/cuda:12.6.2-devel-ubuntu22.04
FROM ${CUDA_IMAGE}

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV CCACHE_DIR=/home/devuser/.ccache
ENV CCACHE_MAXSIZE=1G
ENV CXX_STANDARD=17

# Install essential packages (nvcc + host gcc ship in the CUDA devel image).
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    wget \
    curl \
    fish \
    zsh \
    ninja-build \
    ccache \
    rsync \
    python3 \
    python3-pip \
    python3-dev \
    libgmp-dev \
    libssl-dev \
    sudo \
    software-properties-common \
    lsb-release \
    gnupg \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Create a new user with sudo privileges
RUN useradd -m devuser \
    && echo "devuser:password" | chpasswd \
    && usermod -aG sudo devuser \
    && echo "devuser ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers \
    && mkdir -p $CCACHE_DIR \
    && chown devuser:devuser $CCACHE_DIR

# Set up Python tools
RUN pip3 install --no-cache-dir --upgrade pip setuptools wheel pre-commit

# Add Kitware APT repository for a recent CMake (project requires >= 3.24)
RUN wget -qO- https://apt.kitware.com/keys/kitware-archive-latest.asc | \
    gpg --dearmor -o /usr/share/keyrings/kitware-archive-keyring.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main" | \
    tee /etc/apt/sources.list.d/kitware.list > /dev/null \
    && apt-get update && apt-get install -y --no-install-recommends cmake \
    && rm -rf /var/lib/apt/lists/* \
    && cmake --version && nvcc --version

# Install LLVM/Clang and clang-format 18 (matches the default dev container so
# pre-commit's clang-format hook behaves identically).
RUN wget -q https://apt.llvm.org/llvm.sh -O /tmp/llvm.sh \
    && chmod +x /tmp/llvm.sh \
    && /tmp/llvm.sh 18 || true
RUN apt-get update \
    && apt-get install -y --no-install-recommends clang-18 clang-tools-18 clang-format-18 \
    && rm -rf /var/lib/apt/lists/* \
    && clang-format-18 --version

# Set the default user and working directory
USER devuser
WORKDIR /home/devuser/workspace

CMD ["bash"]
