cmake_minimum_required(VERSION 3.14)

# Project name
project(mirco)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wall -Wextra -O3 -std=c++17")

# Add libraries
add_subdirectory("extern/googletest")

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
include_directories(SYSTEM ${Trilinos_TPL_INCLUDE_DIRS})
include_directories(SYSTEM ${Trilinos_INCLUDE_DIRS})

# Add to the library path Trilinos' library path, and the library
# paths of the third-party libraries (TPLs) with which Trilinos was
# built.
link_directories(${Trilinos_TPL_LIBRARY_DIRS})
link_directories(${Trilinos_LIBRARY_DIRS})

#find JsonCpp package
find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)
link_libraries(${JSONCPP_LIBRARIES})

# Compile mirco library
add_library(mirco_library
  src/topology.cpp
  src/linearsolver.cpp
  src/nonlinearsolver.cpp
  src/matrixsetup.cpp
  src/warmstart.cpp
  src/evaluate.cpp
  src/topologyfactory.cpp
  src/surfmaxandmean.cpp
  src/secondpredictor.cpp
  src/computeresidual.cpp
  src/firstpredictor.cpp
  src/createmeshgrid.cpp
  src/computecontactnodes.cpp
  src/contactforceandarea.cpp
)

add_library(mirco_io
  src/setparameters.cpp
  src/writetofile.cpp
)

add_library(mirco_utils
  src/filesystem_utils.cpp
)

#find OpenMP
find_package(OpenMP REQUIRED)

#find OpenMPI
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

# Compile mirco itself
add_executable(mirco src/main.cpp)

target_link_libraries(mirco
  PUBLIC OpenMP::OpenMP_CXX
  mirco_library
  mirco_io
  mirco_utils
  ${JSONCPP_LIBRARIES}
  ${Trilinos_LIBRARIES}
  ${Trilinos_TPL_LIBRARIES}
  ${MPI_LIBRARIES}
  )

enable_testing()
include(GoogleTest)

#Compile unittest executable
add_executable(tests
  tests/unittests/test.cpp
  tests/unittests/nonlinear_solver_test.cpp
)

target_link_libraries(tests
  PUBLIC OpenMP::OpenMP_CXX
  mirco_library
  mirco_io
  mirco_utils
  gtest
  ${JSONCPP_LIBRARIES}
  ${Trilinos_LIBRARIES}
  ${Trilinos_TPL_LIBRARIES}
  ${MPI_LIBRARIES}
  )

include(GoogleTest)
gtest_discover_tests(tests)
