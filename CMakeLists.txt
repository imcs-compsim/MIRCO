cmake_minimum_required(VERSION 3.14)

# Project name
project(mirco VERSION 0.1.0)

# Check for out-of-source build
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "\nIn-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n\nTo clean-up, you may need to remove CMakeCache.txt and CMakeFiles from your source directory.")
endif()

set(CMAKE_CXX_STANDARD 17)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3")

# Add libraries
if(GTEST_IN_MIRCO)
  add_subdirectory("extern/googletest")
endif()

option(TRILINOS_IN_MIRCO "Find the Trilinos library in the MIRCO project. If set to OFF, use a library from an upstream project" ON)

if(TRILINOS_IN_MIRCO)
  set(CMAKE_PREFIX_PATH ${TRILINOS_PATH} ${CMAKE_PREFIX_PATH})
  find_package(Trilinos REQUIRED HINTS ${TRILINOS_PATH} COMPONENTS Teuchos)

  # handle the case of Trililnos < 14.0 by faking the interface target
  if(Trilinos_FOUND AND NOT TARGET Trilinos::all_selected_libs)
    # In preparation for newer Trilinos releases, create a target
    # Trilinos::all_selected_libs with the correct dependencies
    add_library(_Trilinos_all_selected_libs INTERFACE)
    target_include_directories(
      _Trilinos_all_selected_libs
      SYSTEM
      INTERFACE ${Trilinos_INCLUDE_DIRS}
      INTERFACE ${Trilinos_TPL_INCLUDE_DIRS}
    )
    target_link_libraries(
      _Trilinos_all_selected_libs INTERFACE ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
    )
    add_library(Trilinos::all_selected_libs ALIAS _Trilinos_all_selected_libs)
  endif()
else()
  # Use a Trilinos which has been found in an upstream project. This mechanism only supports the interface target.
  if(TARGET Trilinos::all_selected_libs)
    message(STATUS "Trilinos included from upstream project.")
    add_library(_Trilinos_all_selected_libs INTERFACE)
    target_link_libraries(
      _Trilinos_all_selected_libs INTERFACE Trilinos::all_selected_libs
    )
  else()
    message(FATAL_ERROR "TRILINOS_IN_MIRCO=OFF but no Trilinos target was available in the upstream project.")
  endif()
endif()

# Print details on the Trilinos configuration
if(TRILINOS_IN_MIRCO)
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
endif()

# Find OpenMP
find_package(OpenMP REQUIRED)

# Find OpenMPI
find_package(MPI REQUIRED)

# Compile mirco library
add_library(mirco_core
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
target_link_libraries(mirco_core PUBLIC Trilinos::all_selected_libs OpenMP::OpenMP_CXX)

add_library(mirco_io
  src/mirco_setparameters.cpp
  )
target_link_libraries(mirco_io PUBLIC Trilinos::all_selected_libs)

add_library(mirco_utils
  src/mirco_filesystem_utils.cpp
  )

add_library(mirco_lib INTERFACE)
target_link_libraries(mirco_lib INTERFACE mirco_core mirco_io mirco_utils)
target_link_libraries(mirco_lib
  INTERFACE OpenMP::OpenMP_CXX
  Trilinos::all_selected_libs
  MPI::MPI_CXX
  )

target_include_directories(mirco_lib INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/>
  $<INSTALL_INTERFACE:include/mirco/>
  )

add_library(mirco::mirco_lib ALIAS mirco_lib)

# Compile mirco itself
add_executable(mirco src/main.cpp)
target_link_libraries(mirco PUBLIC mirco::mirco_lib)

# Install mirco (to be used as a library by other codes)
install(TARGETS mirco_lib mirco_core mirco_io mirco_utils _Trilinos_all_selected_libs
  EXPORT mirco_libTargets
  ARCHIVE LIBRARY PUBLIC_HEADER
  )

install(EXPORT mirco_libTargets
  NAMESPACE mirco::
  DESTINATION lib/cmake/mirco
  )

include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/mirco_libConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/mirco_libConfig.cmake
  INSTALL_DESTINATION lib/cmake/mirco
  )

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/mirco_libConfig.cmake
  DESTINATION lib/cmake/mirco 
  )

install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/ DESTINATION include/mirco FILES_MATCHING PATTERN "*.h*")

option(GTEST_IN_MIRCO "Enable GoogleTest in MIRCO" ON)
# Compile unittest executable
if(GTEST_IN_MIRCO)
  enable_testing()
  include(GoogleTest)

  add_executable(tests
    tests/unittests/test.cpp
    tests/unittests/nonlinear_solver_test.cpp
    )

  target_link_libraries(tests PUBLIC mirco::mirco_lib gtest)

  include(GoogleTest)
  gtest_discover_tests(tests)

  include(TestingFramework.cmake)
endif()
