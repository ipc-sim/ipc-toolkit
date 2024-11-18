FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV CCACHE_DIR=/home/devuser/.ccache
ENV CCACHE_MAXSIZE=1G
ENV CXX_STANDARD=17

# Update and install basic dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    wget \
    curl \
    fish \
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
    fish \
    zsh \
    && rm -rf /var/lib/apt/lists/*

# Create a new user with sudo privileges
RUN useradd -m devuser && echo "devuser:password" | chpasswd && usermod -aG sudo devuser && \
    echo "devuser ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers && \
    mkdir -p $CCACHE_DIR && chown devuser:devuser $CCACHE_DIR

# Set up Python tools
RUN pip3 install --upgrade pip setuptools wheel pre-commit

# Install CMake from Kitware APT repository
RUN wget -qO- https://apt.kitware.com/keys/kitware-archive-latest.asc | \
    gpg --dearmor -o /usr/share/keyrings/kitware-archive-keyring.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main" | \
    tee /etc/apt/sources.list.d/kitware.list > /dev/null && \
    apt-get update && \
    apt-get install -y cmake && \
    cmake --version

# Install LLVM/Clang and Clang-Format version 18
RUN wget -q https://apt.llvm.org/llvm.sh -O /tmp/llvm.sh && \
    chmod +x /tmp/llvm.sh && \
    /tmp/llvm.sh 18 || ( \
        echo "llvm.sh script failed, attempting manual installation." && \
        apt-get update && \
        apt-get install -y clang-18 clang-tools-18 clang-format-18 || { \
            echo "Manual LLVM installation failed." && exit 1; \
        } \
    ) && \
    ls -l /usr/bin/clang-format* && ls -l /usr/bin/clang* && \
    ln -s /usr/bin/clang-format-18 /usr/bin/clang-format && \
    ln -s /usr/bin/clang++-18 /usr/bin/clang++ && \
    clang-18 --version && clang-format --version && clang++ --version

# Cleanup to reduce image size
RUN apt-get remove -y git wget curl && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /root/src

# Set the default user and working directory
USER devuser
WORKDIR /home/devuser/workspace

CMD ["bash"]
