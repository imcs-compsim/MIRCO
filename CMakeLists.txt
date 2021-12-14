cmake_minimum_required(VERSION 3.14)

# Project name
project(bem)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -std=c++17")

# Add libraries
set(CMAKE_PREFIX_PATH ${TRILINOS_PATH} ${CMAKE_PREFIX_PATH})
find_package(Trilinos REQUIRED)

# Print details on the Trilinos configuration
if(Trilinos_FOUND)
   message("\nFound Trilinos!  Here are the details: ")
   message("   Trilinos_DIR = ${Trilinos_DIR}")
   message("   Trilinos_VERSION = ${Trilinos_VERSION}")
   message("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
   message("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
   message("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
   message("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
   message("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
   message("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
   message("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
   message("   Trilinos_CXX_COMPILER = ${Trilinos_CXX_COMPILER}")
   message("   Trilinos_C_COMPILER = ${Trilinos_C_COMPILER}")
   message("   Trilinos_Fortran_COMPILER = ${Trilinos_Fortran_COMPILER}")
   message("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
   message("   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
   message("   Trilinos_Fortran_COMPILER_FLAGS = ${Trilinos_Fortran_COMPILER_FLAGS}")
   message("   Trilinos_LINKER = ${Trilinos_LINKER}")
   message("   Trilinos_EXTRA_LD_FLAGS = ${Trilinos_EXTRA_LD_FLAGS}")
   message("   Trilinos_AR = ${Trilinos_AR}")
   message("End of Trilinos details\n")
else()
  message(FATAL_ERROR "Could not find Trilinos!")
endif()

# Specify include directories in this repository and for Trilinos
include_directories(".")
include_directories(${Trilinos_TPL_INCLUDE_DIRS})
include_directories(${Trilinos_INCLUDE_DIRS})

# Add to the library path Trilinos' library path, and the library
# paths of the third-party libraries (TPLs) with which Trilinos was
# built.
link_directories(${Trilinos_TPL_LIBRARY_DIRS})
link_directories(${Trilinos_LIBRARY_DIRS})

# Compile bem code
add_executable(bem src/main.cpp)

#find JsonCpp package
find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)
link_libraries(${JSONCPP_LIBRARIES})

# Compile bem library
add_library(bem_library
  src/topology.cpp
  src/linearsolver.cpp
  src/nonlinearsolver.cpp
  src/matrixsetup.cpp
  src/warmstart.cpp
  src/setparameters.cpp
  src/evaluate.cpp
  src/writetofile.cpp
  src/topologyfactory.cpp
)

#find OpenMP
find_package(OpenMP REQUIRED)

#find OpenMPI
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

target_link_libraries(bem
  PUBLIC OpenMP::OpenMP_CXX
  bem_library
  ${JSONCPP_LIBRARIES}
  ${Trilinos_LIBRARIES}
  ${Trilinos_TPL_LIBRARIES}
  ${MPI_LIBRARIES}
  )
