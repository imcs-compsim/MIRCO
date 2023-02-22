cmake_minimum_required(VERSION 3.14)

# Project name
project(mirco VERSION 0.1.0)

# Check for out-of-source build
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "\nIn-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n\nTo clean-up, you may need to remove CMakeCache.txt and CMakeFiles from your source directory.")
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wall -Wextra -Werror -O3 -std=c++17")

# Add libraries
if(GTEST_IN_MIRCO)
  add_subdirectory("extern/googletest")
endif()

if(TRILINOS_IN_MIRCO)
  set(CMAKE_PREFIX_PATH ${TRILINOS_PATH} ${CMAKE_PREFIX_PATH})
  find_package(Trilinos REQUIRED COMPONENTS Teuchos)
endif()

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

# Compile mirco library
add_library(mirco_library
  src/mirco_topology.cpp
  src/mirco_linearsolver.cpp
  src/mirco_nonlinearsolver.cpp
  src/mirco_matrixsetup.cpp
  src/mirco_warmstart.cpp
  src/mirco_evaluate.cpp
  src/mirco_topologyutilities.cpp
  src/mirco_contactpredictors.cpp
  src/mirco_contactstatus.cpp
)

add_library(mirco_io
  src/mirco_setparameters.cpp
)

add_library(mirco_utils
  src/mirco_filesystem_utils.cpp
)

# Find OpenMP
find_package(OpenMP REQUIRED)

# Find OpenMPI
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

# Compile mirco itself
add_executable(mirco src/main.cpp)

# Install mirco (to be used as a library by other codes)
install(TARGETS mirco DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/ DESTINATION include/mirco FILES_MATCHING PATTERN "*.h*")

target_link_libraries(mirco
  PUBLIC OpenMP::OpenMP_CXX
  mirco_library
  mirco_io
  mirco_utils
  ${Trilinos_LIBRARIES}
  ${Trilinos_TPL_LIBRARIES}
  ${MPI_LIBRARIES}
  )

# Compile unittest executable
if(GTEST_IN_MIRCO)
  enable_testing()
  include(GoogleTest)

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
    ${Trilinos_LIBRARIES}
    ${Trilinos_TPL_LIBRARIES}
    ${MPI_LIBRARIES}
  )

  include(GoogleTest)
  gtest_discover_tests(tests)

  include(TestingFramework.cmake)
endif()

