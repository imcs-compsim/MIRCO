#!/bin/bash

# Install Kokkos with OpenMP backend
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

INSTALL_DIR="$1"
# Number of procs for building (default 4)
NPROCS=${NPROCS:=4}
# git sha from Kokkos repository:
VERSION="f11fb0be01b0bcb2d7fa9eda41f2fe79d63e859b"

CMAKE_COMMAND=cmake

git clone https://github.com/kokkos/kokkos.git kokkos
cd kokkos
git checkout $VERSION
cd .. && mkdir kokkos_build && cd kokkos_build

$CMAKE_COMMAND \
  -D CMAKE_BUILD_TYPE:STRING="RELEASE" \
  -D CMAKE_CXX_STANDARD:STRING="17" \
  -D CMAKE_CXX_COMPILER=g++ \
  -D CMAKE_INSTALL_PREFIX:STRING=$INSTALL_DIR \
  -D BUILD_SHARED_LIBS:BOOL=OFF \
  \
  -D Kokkos_ENABLE_OPENMP=ON \
  \
  ../kokkos

make -j${NPROCS} install
cd ..
rm -rf kokkos kokkos_build
