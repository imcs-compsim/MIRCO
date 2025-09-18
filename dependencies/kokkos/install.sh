#!/bin/bash

# Install Kokkos with OpenMP and with Serial+CUDA backends (separately)
# Call with
# ./install.sh /path/to/install/root

# Exit the script at the first failure
set -e

DEPS_ROOT="$1"
# Number of procs for building (default 4)
NPROCS=${NPROCS:=4}
# git sha from Kokkos repository:
VERSION="f11fb0be01b0bcb2d7fa9eda41f2fe79d63e859b"

CMAKE_COMMAND=cmake

cd $DEPS_ROOT
git clone https://github.com/kokkos/kokkos.git kokkos
cd kokkos
git checkout $VERSION
cd .. && mkdir kokkos_build && cd kokkos_build

# Buld and install with OpenMP backend
$CMAKE_COMMAND \
  -D CMAKE_BUILD_TYPE:STRING="RELEASE" \
  -D CMAKE_CXX_STANDARD:STRING="17" \
  -D CMAKE_CXX_COMPILER=g++ \
  -D CMAKE_INSTALL_PREFIX:STRING=$DEPS_ROOT/kokkos_install_openmp \
  -D BUILD_SHARED_LIBS:BOOL=OFF \
  \
  -D Kokkos_ENABLE_OPENMP=ON \
  \
  ../kokkos
make -j${NPROCS} install

# Clean build dir
rm -rf *

# Buld and install with Serial+CUDA backend
$CMAKE_COMMAND \
  -D CMAKE_BUILD_TYPE:STRING="RELEASE" \
  -D CMAKE_CXX_STANDARD:STRING="17" \
  -D CMAKE_CXX_COMPILER=g++ \
  -D CMAKE_INSTALL_PREFIX:STRING=$DEPS_ROOT/kokkos_install_cuda \
  -D BUILD_SHARED_LIBS:BOOL=OFF \
  \
  -D Kokkos_ENABLE_SERIAL=ON \
  -D Kokkos_ENABLE_CUDA=ON \
  \
  ../kokkos
make -j${NPROCS} install

# keep kokkos (src) for nvcc_wrapper
rm -rf kokkos_build
