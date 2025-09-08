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
  src/mirco_nonlinearsolver.cpp
  )
target_link_libraries(mirco_needKK PRIVATE Kokkos::kokkos KokkosKernels::kokkoskernels)

# Compile mirco library
add_library(mirco_core
  src/mirco_evaluate.cpp
  src/mirco_matrixsetup.cpp
  src/mirco_contactpredictors.cpp
  src/mirco_contactstatus.cpp
  src/mirco_warmstart.cpp
  )
target_link_libraries(mirco_core PRIVATE mirco_needKK Kokkos::kokkos)

add_library(mirco_topology
  src/mirco_topology.cpp
  src/mirco_topologyutilities.cpp
  )
target_link_libraries(mirco_topology PRIVATE Kokkos::kokkos)

add_library(mirco_inputparameters
  src/mirco_inputparameters.cpp
  )
target_link_libraries(mirco_inputparameters PRIVATE mirco_topology Kokkos::kokkos)

# Fetch rapidyaml (ryml) directly
#include(FetchContent)
##set(C4CORE_WITH_FASTFLOAT OFF CACHE BOOL "disable fastfloat in c4core" FORCE)
#FetchContent_Declare(rapidyaml GIT_REPOSITORY https://github.com/biojppm/rapidyaml.git GIT_TAG 47ec2fa184209687c20fd5bc05621e1cb1200311)
#FetchContent_MakeAvailable(rapidyaml)










include(ExternalProject)

set(HOST_PREFIX ${CMAKE_BINARY_DIR}/_hostdeps)
set(HOST_CXX /usr/bin/g++-12)

ExternalProject_Add(ryml_build
  GIT_REPOSITORY https://github.com/biojppm/rapidyaml.git
  GIT_TAG 47ec2fa184209687c20fd5bc05621e1cb1200311
  INSTALL_DIR ${HOST_PREFIX}
  CMAKE_ARGS
    -DCMAKE_CXX_COMPILER=${HOST_CXX}
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DCMAKE_POLICY_DEFAULT_CMP0077=NEW
    -DC4CORE_WITH_FASTFLOAT=OFF
    -DC4CORE_INSTALL=ON
    -DRYML_BUILD_TESTS=OFF
    -DRYML_BUILD_BENCHMARKS=OFF
    -DCMAKE_CXX_FLAGS=-w
    -DCMAKE_C_FLAGS=-w
  INSTALL_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target install --config $<CONFIG>
  BUILD_BYPRODUCTS ${HOST_PREFIX}/lib/libryml.a
)

file(MAKE_DIRECTORY ${HOST_PREFIX}/include ${HOST_PREFIX}/lib)

add_library(ryml_import STATIC IMPORTED)
set_target_properties(ryml_import PROPERTIES
  IMPORTED_LOCATION ${HOST_PREFIX}/lib/libryml.a
  INTERFACE_INCLUDE_DIRECTORIES ${HOST_PREFIX}/include
  INTERFACE_COMPILE_DEFINITIONS C4CORE_NO_FAST_FLOAT)
add_dependencies(ryml_import ryml_build)

add_library(ryml_host INTERFACE)
add_dependencies(ryml_host ryml_build)
target_include_directories(ryml_host INTERFACE ${HOST_PREFIX}/include)
target_link_libraries(ryml_host INTERFACE ${HOST_PREFIX}/lib/libryml.a)
target_compile_definitions(ryml_host INTERFACE C4CORE_NO_FAST_FLOAT)





















add_library(mirco_utils
  src/mirco_utils.cpp
  )

target_link_libraries(mirco_utils PUBLIC ryml_host)#ryml)












add_library(mirco_inputparameters_yaml
  src/mirco_inputparameters_yaml.cpp
  )
target_link_libraries(mirco_inputparameters_yaml PRIVATE mirco_utils Kokkos::kokkos)

add_library(mirco_lib INTERFACE)
target_link_libraries(mirco_lib INTERFACE mirco_needKK mirco_core mirco_topology mirco_inputparameters)

# Compile mirco itself
add_executable(mirco src/main.cpp)

target_include_directories(mirco_lib INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/>
  $<INSTALL_INTERFACE:include/mirco/>
  )

add_library(mirco::mirco_lib ALIAS mirco_lib)

target_link_libraries(mirco PUBLIC mirco::mirco_lib mirco_inputparameters_yaml mirco_utils Kokkos::kokkos)

# Install mirco (to be used as a library by other codes)
install(TARGETS mirco_lib mirco_needKK mirco_core mirco_topology mirco_inputparameters
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

install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/ DESTINATION include/mirco FILES_MATCHING PATTERN "*.h")

option(GTEST_IN_MIRCO "Enable GoogleTest in MIRCO" ON)
# Compile unittest executable
if(GTEST_IN_MIRCO)
  enable_testing()
  include(GoogleTest)

  add_executable(tests
    tests/unittests/test.cpp
    )
  target_link_libraries(tests PUBLIC mirco::mirco_lib mirco_inputparameters_yaml mirco_utils gtest Kokkos::kokkos)

  # Copy tests/unittests/data folder to the build dir
  add_custom_command(
    TARGET tests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/unittests/data
            $<TARGET_FILE_DIR:tests>/test/data
  )

  gtest_discover_tests(tests)
  include(TestingFramework.cmake)
endif()
