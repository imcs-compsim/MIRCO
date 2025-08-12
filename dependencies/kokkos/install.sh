#!/bin/bash

# Install trilinos
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

INSTALL_DIR="$1"
# Number of procs for building (default 4)
NPROCS=${NPROCS:=4}
# git sha from Trilinos repository:
VERSION="06db4c850654feacabdaed61ee8308219266b6a5"
#CHECKSUM=""


# Location of script to apply patches later
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
CMAKE_COMMAND=cmake

git clone https://github.com/trilinos/Trilinos.git
cd Trilinos
git checkout $VERSION
cd .. && mkdir trilinos_build && cd trilinos_build

MPI_DIR=/usr
MPI_BIN_DIR=$MPI_DIR/bin

$CMAKE_COMMAND \
  -D CMAKE_BUILD_TYPE:STRING="RELEASE" \
  -D CMAKE_CXX_STANDARD:STRING="17" \
  -D CMAKE_CXX_COMPILER:FILEPATH="$MPI_BIN_DIR/mpic++" \
  -D CMAKE_C_COMPILER:FILEPATH="$MPI_BIN_DIR/mpicc" \
  -D CMAKE_Fortran_COMPILER:FILEPATH="$MPI_BIN_DIR/mpif90" \
  -D CMAKE_INSTALL_PREFIX:STRING=$INSTALL_DIR \
  -D BUILD_SHARED_LIBS:BOOL=ON \
  \
  -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_TESTS:BOOL=OFF \
  -D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
  \
  -D Trilinos_ASSERT_MISSING_PACKAGES=OFF \
  -D Trilinos_ENABLE_Gtest:BOOL=OFF \
  -D Trilinos_ENABLE_Teuchos:BOOL=ON \
  \
  ../Trilinos

make -j${NPROCS} install
cd ..
rm -rf Trilinos trilinos_build
