#!/bin/bash

# Install Kokkos-Kernels with OpenMP backend
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

INSTALL_DIR="$1"
# Number of procs for building (default 4)
NPROCS=${NPROCS:=4}
# git sha from Kokkos-Kernels repository:
VERSION="3254a1c1ccda11673fad64651dd8ab957bf49e7d"

CMAKE_COMMAND=cmake

git clone https://github.com/kokkos/kokkos-kernels.git kokkos-kernels
cd kokkos-kernels
git checkout $VERSION
cd .. && mkdir kokkos-kernels_build && cd kokkos-kernels_build

$CMAKE_COMMAND \
  -D CMAKE_BUILD_TYPE:STRING="RELEASE" \
  -D CMAKE_CXX_STANDARD:STRING="17" \
  -D CMAKE_CXX_COMPILER=g++ \
  -D CMAKE_INSTALL_PREFIX:STRING=$INSTALL_DIR \
  \
  -D Kokkos_ROOT=$INSTALL_DIR/../kokkos_install \
  \
  -D KokkosKernels_ENABLE_TPL_BLAS=ON \
  -D KokkosKernels_ENABLE_TPL_LAPACK=ON \
  \
  -D Kokkos_ENABLE_OPENMP=ON \
  \
  ../kokkos-kernels

make -j${NPROCS} install
cd ..
rm -rf kokkos-kernels kokkos-kernels_build
