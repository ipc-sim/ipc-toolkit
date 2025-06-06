FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV CCACHE_DIR=/home/devuser/.ccache
ENV CCACHE_MAXSIZE=1G
ENV CXX_STANDARD=17

# Update package lists
RUN apt-get update

# Install essential packages
RUN apt-get install -y --no-install-recommends \
    build-essential \
    git \
    wget \
    curl \
    fish \
    zsh \
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
    software-properties-common \
    lsb-release \
    gnupg \
    && rm -rf /var/lib/apt/lists/*

# Create a new user with sudo privileges
RUN useradd -m devuser \
    && echo "devuser:password" | chpasswd \
    && usermod -aG sudo devuser \
    && echo "devuser ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers \
    && mkdir -p $CCACHE_DIR \
    && chown devuser:devuser $CCACHE_DIR

# Set up Python tools
RUN pip3 install --upgrade pip setuptools wheel pre-commit

# Add Kitware APT repository for CMake
RUN wget -qO- https://apt.kitware.com/keys/kitware-archive-latest.asc | \
    gpg --dearmor -o /usr/share/keyrings/kitware-archive-keyring.gpg

RUN echo "deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main" | \
    tee /etc/apt/sources.list.d/kitware.list > /dev/null

# Update package lists and install CMake
RUN apt-get update
RUN apt-get install -y cmake

# Verify CMake installation
RUN cmake --version

# Install LLVM/Clang and Clang-Format version 18
RUN wget -q https://apt.llvm.org/llvm.sh -O /tmp/llvm.sh
RUN chmod +x /tmp/llvm.sh
RUN /tmp/llvm.sh 18 || true

RUN apt-get update && apt-get install -y clang-18 clang-tools-18 clang-format-18
    
RUN clang-18 --version
RUN clang++-18 --version
RUN clang-format-18 --version

# 12. Set the default user and working directory
USER devuser

WORKDIR /home/devuser/workspace

CMD ["bash"]
