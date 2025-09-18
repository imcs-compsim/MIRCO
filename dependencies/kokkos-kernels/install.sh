#!/bin/bash

# Install Kokkos-Kernels with OpenMP and with Serial+CUDA backends (separately)
# Call with
# ./install.sh /path/to/install/root

# Exit the script at the first failure
set -e

INSTALL_ROOT="$1"
# Number of procs for building (default 4)
NPROCS=${NPROCS:=4}
# git sha from Kokkos-Kernels repository:
VERSION="3254a1c1ccda11673fad64651dd8ab957bf49e7d"

CMAKE_COMMAND=cmake

git clone https://github.com/kokkos/kokkos-kernels.git kokkos-kernels
cd kokkos-kernels
git checkout $VERSION
cd .. && mkdir kokkos-kernels_build && cd kokkos-kernels_build

# Buld and install with OpenMP backend
$CMAKE_COMMAND \
  -D CMAKE_BUILD_TYPE:STRING="RELEASE" \
  -D CMAKE_CXX_STANDARD:STRING="17" \
  -D CMAKE_CXX_COMPILER=g++ \
  -D CMAKE_INSTALL_PREFIX:STRING=$INSTALL_ROOT/kokkos-kernels_install_openmp \
  -D BUILD_SHARED_LIBS:BOOL=OFF \
  \
  -D Kokkos_ROOT=$INSTALL_ROOT/../kokkos_install_openmp \
  \
  -D KokkosKernels_ENABLE_TPL_BLAS=ON \
  -D KokkosKernels_ENABLE_TPL_LAPACK=ON \
  \
  -D Kokkos_ENABLE_OPENMP=ON \
  \
  ../kokkos-kernels
make -j${NPROCS} install

# Clean build dir
rm -rf *

# Buld and install with Serial+CUDA backend
$CMAKE_COMMAND \
  -D CMAKE_BUILD_TYPE:STRING="RELEASE" \
  -D CMAKE_CXX_STANDARD:STRING="17" \
  -D CMAKE_CXX_COMPILER=g++ \
  -D CMAKE_INSTALL_PREFIX:STRING=$INSTALL_ROOT/kokkos-kernels_install_cuda \
  -D BUILD_SHARED_LIBS:BOOL=OFF \
  \
  -D Kokkos_ROOT=$INSTALL_ROOT/../kokkos_install_cuda \
  \
  -D KokkosKernels_ENABLE_TPL_BLAS=ON \
  -D KokkosKernels_ENABLE_TPL_LAPACK=ON \
  -D KokkosKernels_ENABLE_TPL_CUSOLVER=ON \
  \
  -D Kokkos_ENABLE_SERIAL=ON \
  -D Kokkos_ENABLE_CUDA=ON \
  \
  ../kokkos-kernels
make -j${NPROCS} install

rm -rf kokkos kokkos_build
