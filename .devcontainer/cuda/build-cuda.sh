#!/usr/bin/env bash
#
# Headless CUDA compile check that reuses the CUDA dev-container image.
#
# The dev container (devcontainer.json) is for *interactive* development. This
# script is the *batch* counterpart: it builds the same image and compiles the
# whole project with CUDA enabled, then exits with the build's status. Use it
# from a Mac (compile-only; no GPU) or in CI to keep the CUDA build green.
#
# The source is mounted read-only and copied into the container (minus build
# artifacts and the machine-specific IPCToolkitOptions.cmake) so the build is
# hermetic and never writes into your host working tree.
#
# Usage:
#   .devcontainer/cuda/build-cuda.sh                 # cuda-release, arch 75;80;86;89
#   PRESET=test .devcontainer/cuda/build-cuda.sh     # test preset (CUDA + tests)
#   CUDA_ARCH="86" .devcontainer/cuda/build-cuda.sh  # single architecture
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

IMAGE_NAME="${IMAGE_NAME:-ipc-toolkit-cuda-dev}"
PRESET="${PRESET:-cuda-release}"
CUDA_ARCH="${CUDA_ARCH:-75;80;86;89}"
CUDA_IMAGE="${CUDA_IMAGE:-nvidia/cuda:12.6.2-devel-ubuntu22.04}"

echo ">> Building CUDA dev image '${IMAGE_NAME}'"
docker build \
    -f "${REPO_ROOT}/.devcontainer/cuda/Dockerfile" \
    -t "${IMAGE_NAME}" \
    --build-arg "CUDA_IMAGE=${CUDA_IMAGE}" \
    "${REPO_ROOT}"

echo ">> Compiling (preset=${PRESET}, arch=${CUDA_ARCH})"
# Run as root so the named cache volumes are writable; the source is mounted
# read-only and copied to a scratch dir inside the container.
docker run --rm --user root \
    -e PRESET="${PRESET}" \
    -e CUDA_ARCH="${CUDA_ARCH}" \
    -v "${REPO_ROOT}":/src:ro \
    -v ipc-toolkit-cpm-cache:/cpm-cache \
    -v ipc-toolkit-cuda-ccache:/root/.ccache \
    "${IMAGE_NAME}" \
    bash -euo pipefail -c '
        export CPM_SOURCE_CACHE=/cpm-cache CCACHE_DIR=/root/.ccache
        mkdir -p /workspace
        tar -C /src --exclude=./build --exclude=./.git -cf - . \
            | tar -C /workspace -xf -
        rm -f /workspace/IPCToolkitOptions.cmake
        cd /workspace
        cmake --preset="${PRESET}" -G Ninja \
            -DCMAKE_CUDA_ARCHITECTURES="${CUDA_ARCH}" \
            -DSCALABLE_CCD_CUDA_ARCHITECTURES="${CUDA_ARCH}" \
            -DCMAKE_CXX_FLAGS="-Wno-psabi" \
            -DCMAKE_CUDA_FLAGS="-Xcompiler=-Wno-psabi"
        cmake --build --preset="${PRESET}" -j "$(nproc)"
    '

echo ">> Done. CUDA build succeeded (compile-only; code was not executed)."
