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
VERSION="b95780f4ea5b50c765046a8b694a4f196837653d"

CMAKE_COMMAND=cmake

cd $DEPS_ROOT
git clone https://github.com/kokkos/kokkos.git kokkos
cd kokkos
git checkout $VERSION
cd .. && mkdir kokkos_build && cd kokkos_build

# Build and install with Serial backend
$CMAKE_COMMAND \
  -D CMAKE_BUILD_TYPE:STRING="RELEASE" \
  -D CMAKE_CXX_STANDARD:STRING="20" \
  -D CMAKE_CXX_COMPILER=g++ \
  -D CMAKE_INSTALL_PREFIX:STRING=$DEPS_ROOT/kokkos_install_serial \
  -D BUILD_SHARED_LIBS:BOOL=OFF \
  \
  -D Kokkos_ENABLE_SERIAL=ON \
  \
../kokkos
make -j${NPROCS} install

# Clean build dir
rm -rf *

# Build and install with OpenMP backend
$CMAKE_COMMAND \
  -D CMAKE_BUILD_TYPE:STRING="RELEASE" \
  -D CMAKE_CXX_STANDARD:STRING="20" \
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

# Build and install with Serial+CUDA backend
$CMAKE_COMMAND \
  -D CMAKE_BUILD_TYPE:STRING="RELEASE" \
  -D CMAKE_CXX_STANDARD:STRING="20" \
  -D CMAKE_CXX_COMPILER=g++ \
  -D CMAKE_INSTALL_PREFIX:STRING=$DEPS_ROOT/kokkos_install_cuda \
  -D BUILD_SHARED_LIBS:BOOL=OFF \
  \
  -D Kokkos_ENABLE_SERIAL=ON \
  -D Kokkos_ENABLE_CUDA=ON \
  \
  -D Kokkos_ARCH_HOPPER90=ON \
  \
../kokkos
make -j${NPROCS} install

# Clean build dir
rm -rf *

# keep kokkos (src) for nvcc_wrapper
cd .. && rm -rf kokkos_build
