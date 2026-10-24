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
#   .devcontainer/cuda/build-cuda.sh                 # cuda-release, arch 75
#   PRESET=test .devcontainer/cuda/build-cuda.sh     # test preset (CUDA + tests)
#   CUDA_ARCH="75;80;86;89" .devcontainer/cuda/build-cuda.sh  # several archs
#   JOBS=4 .devcontainer/cuda/build-cuda.sh          # limit parallelism (memory)
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

IMAGE_NAME="${IMAGE_NAME:-ipc-toolkit-cuda-dev}"
PRESET="${PRESET:-cuda-release}"
# One architecture is enough to answer "does it compile": nvcc runs the whole
# device front-end per architecture in the list, so the extra ones only repeat
# codegen. 75 is the oldest we support, hence the strictest. Override to build
# a list when you want to check architecture-specific codegen.
CUDA_ARCH="${CUDA_ARCH:-75}"
CUDA_IMAGE="${CUDA_IMAGE:-nvidia/cuda:12.6.2-devel-ubuntu22.04}"
# Heavy TUs (headers textually include implementations under CUDA) can OOM the
# VM at full parallelism; default below nproc.
JOBS="${JOBS:-4}"

echo ">> Building CUDA dev image '${IMAGE_NAME}'"
docker build \
    -f "${REPO_ROOT}/.devcontainer/Dockerfile" \
    -t "${IMAGE_NAME}" \
    --build-arg "BASE_IMAGE=${CUDA_IMAGE}" \
    "${REPO_ROOT}"

echo ">> Compiling (preset=${PRESET}, arch=${CUDA_ARCH})"
# Run as root so the named cache volumes are writable; the source is mounted
# read-only and copied to a scratch dir inside the container.
docker run --rm --user root \
    -e PRESET="${PRESET}" \
    -e CUDA_ARCH="${CUDA_ARCH}" \
    -e JOBS="${JOBS}" \
    -v "${REPO_ROOT}":/src:ro \
    -v ipc-toolkit-cuda-workspace:/workspace \
    -v ipc-toolkit-cpm-cache:/cpm-cache \
    -v ipc-toolkit-cuda-ccache:/root/.ccache \
    "${IMAGE_NAME}" \
    bash -euo pipefail -c '
        export CPM_SOURCE_CACHE=/cpm-cache CCACHE_DIR=/root/.ccache
        # /workspace is a persistent named volume: rsync copies only files
        # that changed since the last run (the macOS<->VM file-share is slow,
        # so minimizing reads matters) and ninja can then build incrementally.
        # The excludes also shield the persistent build/ dir from --delete.
        # tests/data is excluded for the same reason: it is cloned by an
        # ExternalProject whose stamp lives under build/, so deleting the data
        # while keeping the stamp makes the next build fail in gitupdate.cmake
        # instead of re-cloning.
        echo ">> [1/3] Syncing source into the container (delta copy)..."
        time rsync -a --delete \
            --exclude=/build \
            --exclude=/.git \
            --exclude=/.ccache \
            --exclude=/tests/data \
            --exclude=/docs \
            --exclude=/notebooks \
            --exclude=/IPCToolkitOptions.cmake \
            /src/ /workspace/
        cd /workspace
        echo ">> [2/3] Configuring (preset=${PRESET})..."
        cmake --preset="${PRESET}" -G Ninja \
            -DCMAKE_CUDA_ARCHITECTURES="${CUDA_ARCH}" \
            -DSCALABLE_CCD_CUDA_ARCHITECTURES="${CUDA_ARCH}" \
            -DCMAKE_CXX_FLAGS="-Wno-psabi" \
            -DCMAKE_CUDA_FLAGS="-Xcompiler=-Wno-psabi"
        echo ">> [3/3] Building (-j ${JOBS})..."
        cmake --build --preset="${PRESET}" -j "${JOBS}"
    '

echo ">> Done. CUDA build succeeded (compile-only; code was not executed)."
