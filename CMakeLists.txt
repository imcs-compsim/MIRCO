cmake_minimum_required(VERSION 3.28)

# try to prevent modification of source directory
# note: some files may still be written before CMake can abort and need to be removed manually
if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  message(
    FATAL_ERROR
      "In-source build not allowed. "
      "Please create a new directory, preferably next to the source directory, and run CMake from there. "
      "You may want to remove CMakeCache.txt and CMakeFiles/ which were created in the source directory."
    )
endif()

project(mirco VERSION 0.1.0)

# Print CMake version to screen
message(STATUS "Using CMake ${CMAKE_VERSION}")

# Enforce the C++ standard we are using and turn off compiler-specific extensions
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# We do not use C++ modules (yet). Turn off scanning to avoid issues with clang-tidy.
# If you want to add module support, this problem needs to be revisited. Our hope
# is that CMake and/or clang-tidy will be updated to handle modules better.
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

# Ensure cmake setup the correct runtime path when installing
# see here for more information: https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling
include(GNUInstallDirs)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Add libraries
if(GTEST_IN_MIRCO)
  add_subdirectory("extern/googletest")
endif()

option(USE_KOKKOS_CODE "Use the Kokkos variant of MIRCO, which also requires Kokkos-Kernels. If set to OFF, the standard OpenMP-based implementation will be used" ON)

if(NOT USE_KOKKOS_CODE)

  option(TRILINOS_IN_MIRCO "Find the Trilinos library in the MIRCO project. If set to OFF, use a library from an upstream project" ON)

  if(TRILINOS_IN_MIRCO)
    set(CMAKE_PREFIX_PATH ${TRILINOS_PATH} ${CMAKE_PREFIX_PATH})
    find_package(Trilinos REQUIRED HINTS ${TRILINOS_PATH} COMPONENTS Teuchos)

    # handle the case of Trililnos < 14.0 by faking the interface target
    if(Trilinos_FOUND AND NOT TARGET Trilinos::all_selected_libs)
      # In preparation for newer Trilinos releases, create a target
      # Trilinos::all_selected_libs with the correct dependencies
      add_library(Trilinos::all_selected_libs IMPORTED INTERFACE)
      target_include_directories(
        Trilinos::all_selected_libs
        SYSTEM
        INTERFACE ${Trilinos_INCLUDE_DIRS}
        INTERFACE ${Trilinos_TPL_INCLUDE_DIRS}
      )
      target_link_libraries(
        Trilinos::all_selected_libs INTERFACE ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
      )
    endif()
  else()
    # Use a Trilinos which has been found in an upstream project. This mechanism only supports the interface target.
    if(TARGET Trilinos::all_selected_libs)
      message(STATUS "Trilinos included from upstream project.")
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

  # Compiler flags: delay this modification after MPI was found. Otherwise,
  # trying to compile with MPI can give warnings.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3")

  # Compile mirco library
  add_library(mirco_core
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

  add_library(mirco_topology
    src/mirco_topology.cpp
    )
  target_link_libraries(mirco_topology PRIVATE Trilinos::all_selected_libs)

  add_library(mirco_utils
    src/mirco_filesystem_utils.cpp
    )

  add_library(mirco_io
    src/mirco_inputparameters.cpp
    )
  target_link_libraries(mirco_io PRIVATE mirco_utils mirco_topology)
  target_link_libraries(mirco_io PUBLIC Trilinos::all_selected_libs)

  add_library(mirco_lib INTERFACE)
  target_link_libraries(mirco_lib INTERFACE mirco_core mirco_io mirco_utils mirco_topology)
  target_link_libraries(mirco_lib
    INTERFACE OpenMP::OpenMP_CXX
    Trilinos::all_selected_libs
    MPI::MPI_CXX
    )

  # All standard/original header files (i.e. without _kokkos)
  file(GLOB_RECURSE HEADER_FILES
    "${PROJECT_SOURCE_DIR}/src/*.h"
  )
  foreach(f ${HEADER_FILES})
    if(${f} MATCHES ".*_kokkos.*")
      list(REMOVE_ITEM HEADER_FILES ${f})
    endif()
  endforeach()

  # Compile mirco itself
  add_executable(mirco src/main.cpp)

else() #if(USE_KOKKOS_CODE)
  # Kokkos
  option(KOKKOS_IN_MIRCO "Find Kokkos in the MIRCO project. If set to OFF, use the Kokkos installation from an upstream project" ON)
  if(KOKKOS_IN_MIRCO)
    set(CMAKE_PREFIX_PATH ${KOKKOS_PATH} ${CMAKE_PREFIX_PATH})
    find_package(Kokkos REQUIRED)
  else()
    # Use a Kokkos from an upstream project
    if(TARGET Kokkos::kokkos)
      message(STATUS "Kokkos included from upstream project.")
    else()
      message(FATAL_ERROR "KOKKOS_IN_MIRCO=OFF but no Kokkos target was available in the upstream project.")
    endif()
  endif()

  # Print details on the Kokkos configuration
  if(KOKKOS_IN_MIRCO)
    if(Kokkos_FOUND)
      message("\nFound Kokkos! Here are the details:")
      get_target_property(Kokkos_INCLUDE_DIRS Kokkos::kokkos INTERFACE_INCLUDE_DIRECTORIES)
      message("   Kokkos_INCLUDE_DIRS = ${Kokkos_INCLUDE_DIRS}")
      get_target_property(Kokkos_COMPILE_DEFINITIONS Kokkos::kokkos INTERFACE_COMPILE_DEFINITIONS)
      message("   Kokkos_COMPILE_DEFINITIONS = ${Kokkos_COMPILE_DEFINITIONS}")
      get_target_property(Kokkos_COMPILE_OPTIONS Kokkos::kokkos INTERFACE_COMPILE_OPTIONS)
      message("   Kokkos_COMPILE_OPTIONS = ${Kokkos_COMPILE_OPTIONS}")
      message("End of Kokkos details\n")
    else()
      message(FATAL_ERROR "Could not find Kokkos!")
    endif()
  endif()

  # Kokkos-Kernels
  option(KOKKOS_KERNELS_IN_MIRCO "Find Kokkos-Kernels in the MIRCO project. If set to OFF, use the Kokkos-Kernels installation from an upstream project" ON)
  if(KOKKOS_KERNELS_IN_MIRCO)
    set(CMAKE_PREFIX_PATH ${KOKKOS_KERNELS_PATH} ${CMAKE_PREFIX_PATH})
    find_package(KokkosKernels REQUIRED)
  else()
  # Use a KokkosKernels from an upstream project
    if(TARGET KokkosKernels::kokkoskernels)
      message(STATUS "KokkosKernels included from upstream project.")
    else()
      message(FATAL_ERROR "KOKKOS_KERNELS_IN_MIRCO=OFF but no KokkosKernels target was available in the upstream project.")
    endif()
  endif()

  # Print details on the Kokkos-Kernels configuration
  if(KOKKOS_KERNELS_IN_MIRCO)
    if(KokkosKernels_FOUND)
      message("\nFound KokkosKernels! Here are the details:")
      get_target_property(KokkosKernels_INCLUDE_DIRS KokkosKernels::kokkoskernels INTERFACE_INCLUDE_DIRECTORIES)
      message("   KokkosKernels_INCLUDE_DIRS = ${KokkosKernels_INCLUDE_DIRS}")
      get_target_property(KokkosKernels_COMPILE_DEFINITIONS KokkosKernels::kokkoskernels INTERFACE_COMPILE_DEFINITIONS)
      message("   KokkosKernels_COMPILE_DEFINITIONS = ${KokkosKernels_COMPILE_DEFINITIONS}")
      get_target_property(KokkosKernels_COMPILE_OPTIONS KokkosKernels::kokkoskernels INTERFACE_COMPILE_OPTIONS)
      message("   KokkosKernels_COMPILE_OPTIONS = ${KokkosKernels_COMPILE_OPTIONS}")
      message("End of KokkosKernels details\n")
    else()
      message(FATAL_ERROR "Could not find KokkosKernels!")
    endif()
  endif()

  # File(s) which need(s) Kokkos-Kernels
  add_library(mirco_needKK
    src/mirco_nonlinearsolver_kokkos.cpp
    )
  target_link_libraries(mirco_needKK PRIVATE Kokkos::kokkos KokkosKernels::kokkoskernels)

  # Compile mirco library
  add_library(mirco_core
    src/mirco_evaluate_kokkos.cpp
    src/mirco_matrixsetup_kokkos.cpp
    src/mirco_contactpredictors_kokkos.cpp
    src/mirco_contactstatus_kokkos.cpp
    src/mirco_warmstart_kokkos.cpp
    )
  target_link_libraries(mirco_core PUBLIC Kokkos::kokkos mirco_needKK)

  add_library(mirco_topology
    src/mirco_topology_kokkos.cpp
    src/mirco_topologyutilities_kokkos.cpp
    )
  target_link_libraries(mirco_topology PUBLIC Kokkos::kokkos)

  add_library(mirco_inputparameters
    src/mirco_inputparameters_kokkos.cpp
    )
  target_link_libraries(mirco_inputparameters PRIVATE mirco_topology)

  # Fetch rapidyaml (ryml) directly
  include(FetchContent)
  FetchContent_Declare(rapidyaml GIT_REPOSITORY https://github.com/biojppm/rapidyaml.git GIT_TAG 47ec2fa184209687c20fd5bc05621e1cb1200311)
  FetchContent_MakeAvailable(rapidyaml)

  add_library(mirco_utils
    src/mirco_utils_kokkos.cpp
    )
  target_link_libraries(mirco_utils PUBLIC ryml)

  add_library(mirco_inputparameters_yaml
    src/mirco_inputparameters_yaml_kokkos.cpp
    )
  target_link_libraries(mirco_inputparameters_yaml PRIVATE mirco_utils Kokkos::kokkos)

  add_library(mirco_lib INTERFACE)
  target_link_libraries(mirco_lib INTERFACE mirco_needKK mirco_core mirco_topology mirco_inputparameters)

  # All _kokkos header files
  file(GLOB_RECURSE HEADER_FILES
    "${PROJECT_SOURCE_DIR}/src/*_kokkos.h"
  )

  # Compile mirco itself
  add_executable(mirco src/main_kokkos.cpp)

endif()

target_include_directories(mirco_lib INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/>
  $<INSTALL_INTERFACE:include/mirco/>
  )

add_library(mirco::mirco_lib ALIAS mirco_lib)

if(NOT USE_KOKKOS_CODE)
  target_link_libraries(mirco PUBLIC mirco::mirco_lib)

  # Install mirco (to be used as a library by other codes)
  install(TARGETS mirco_lib mirco_core mirco_io mirco_utils mirco_topology
    EXPORT mirco_libTargets
    ARCHIVE LIBRARY PUBLIC_HEADER
    )
else()
  target_link_libraries(mirco PUBLIC mirco::mirco_lib mirco_inputparameters_yaml mirco_utils)

  # Install mirco (to be used as a library by other codes)
  install(TARGETS mirco_lib mirco_needKK mirco_core mirco_topology mirco_inputparameters
    EXPORT mirco_libTargets
    ARCHIVE LIBRARY PUBLIC_HEADER
    )
endif()

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

install(FILES ${HEADER_FILES} DESTINATION include/mirco)

option(GTEST_IN_MIRCO "Enable GoogleTest in MIRCO" ON)
# Compile unittest executable
if(GTEST_IN_MIRCO)
  enable_testing()
  include(GoogleTest)

  if(NOT USE_KOKKOS_CODE)
    add_executable(tests
      tests/unittests/test.cpp
      tests/unittests/nonlinear_solver_test.cpp
      )
    target_link_libraries(tests PUBLIC mirco::mirco_lib gtest)
  else()
    add_executable(tests
      tests/unittests/test_kokkos.cpp
      )
    target_link_libraries(tests PUBLIC mirco::mirco_lib mirco_inputparameters_yaml mirco_utils gtest)
  endif()

  # Copy tests/unittests/data folder to the build dir
  add_custom_command(
    TARGET tests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/unittests/data
            $<TARGET_FILE_DIR:tests>/test/data
  )

  gtest_discover_tests(tests)

  if(NOT USE_KOKKOS_CODE)
    include(TestingFramework.cmake)
  else()
    include(TestingFramework_kokkos.cmake)
  endif()
endif()
