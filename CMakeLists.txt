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
set(CMAKE_CXX_STANDARD 20)
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

option(MIRCO_ENABLE_FLAT "Build the flatMirco utility, which can be used to compute shape factors" OFF)

# Add libraries
if(GTEST_IN_MIRCO)
  add_subdirectory("extern/googletest")
endif()

# Kokkos
option(KOKKOS_IN_MIRCO "Find Kokkos in the MIRCO project. If set to OFF, use the Kokkos installation from an upstream project" OFF)
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
option(KOKKOS_KERNELS_IN_MIRCO "Find Kokkos-Kernels in the MIRCO project. If set to OFF, use the Kokkos-Kernels installation from an upstream project" OFF)
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

option(MIRCO_ENABLE_VISUALIZATIONEXPORT "Enable exporting the results for visualization. This requires VTK, either explicitly provided through VTK_PATH, or from an upstream project." OFF)
if(MIRCO_ENABLE_VISUALIZATIONEXPORT)
  # If user provided a prefix, help CMake find VTK there
  if(DEFINED VTK_PATH)
    set(CMAKE_PREFIX_PATH ${VTK_PATH} ${CMAKE_PREFIX_PATH})
  endif()

  find_package(VTK QUIET COMPONENTS CommonCore CommonDataModel IOCore IOXML REQUIRED)

  if(NOT VTK_FOUND)
    message(WARNING "MIRCO_ENABLE_VISUALIZATIONEXPORT is ON, but VTK could not be found!")
    set(MIRCO_WITH_VTK OFF CACHE BOOL "" FORCE)
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

if(MIRCO_ENABLE_VISUALIZATIONEXPORT)
  target_sources(mirco_core PRIVATE src/mirco_exportvisualization.cpp)
  target_link_libraries(mirco_core PRIVATE
    VTK::CommonCore
    VTK::CommonDataModel
    VTK::IOCore
    VTK::IOXML
  )
endif()

add_library(mirco_topology
  src/mirco_topology.cpp
  src/mirco_topologyutilities.cpp
  )
target_link_libraries(mirco_topology PRIVATE Kokkos::kokkos)

add_library(mirco_shapefactors
  src/mirco_shapefactors.cpp
  )

add_library(mirco_inputparameters
  src/mirco_inputparameters.cpp
  )
target_link_libraries(mirco_inputparameters PRIVATE mirco_topology mirco_shapefactors Kokkos::kokkos)

option(RYML_IN_MIRCO "Find ryml in MIRCO project. If set to OFF, use the ryml installation from an upstream project" OFF)

if(RYML_IN_MIRCO)
  message(STATUS "Building ryml within MIRCO project.")
  if(Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP OR Kokkos_ENABLE_SYCL OR Kokkos_ENABLE_OPENMPTARGET)
    # Compile rapidyaml (ryml) separately to avoid conflicts with device-side/offloaded compilers (e.g. nvcc_wrapper)
    message(STATUS "Building rapidyaml using MIRCO_HOST_CXX_COMPILER")
    include(ExternalProject)
    set(RYML_INSTALL_DIR ${CMAKE_BINARY_DIR}/extern/ryml)

    if(NOT DEFINED MIRCO_HOST_CXX_COMPILER)
      set(MIRCO_HOST_CXX_COMPILER "g++-12")
    endif()

    ExternalProject_Add(ryml_build
      GIT_REPOSITORY https://github.com/biojppm/rapidyaml.git
      GIT_TAG 47ec2fa184209687c20fd5bc05621e1cb1200311
      INSTALL_DIR ${RYML_INSTALL_DIR}
      CMAKE_ARGS
        -DCMAKE_CXX_COMPILER=${MIRCO_HOST_CXX_COMPILER}
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DCMAKE_POLICY_DEFAULT_CMP0077=NEW
        -DC4CORE_WITH_FASTFLOAT=OFF
        -DC4CORE_INSTALL=ON
        -DRYML_BUILD_TESTS=OFF
        -DRYML_BUILD_BENCHMARKS=OFF
        -DCMAKE_CXX_FLAGS=-w
      INSTALL_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target install --config $<CONFIG>
      BUILD_BYPRODUCTS ${RYML_INSTALL_DIR}/lib/libryml.a
    )

    file(MAKE_DIRECTORY ${RYML_INSTALL_DIR}/include ${RYML_INSTALL_DIR}/lib)

    add_library(ryml_import STATIC IMPORTED)
    set_target_properties(ryml_import PROPERTIES
      IMPORTED_LOCATION ${RYML_INSTALL_DIR}/lib/libryml.a
      INTERFACE_INCLUDE_DIRECTORIES ${RYML_INSTALL_DIR}/include
      INTERFACE_COMPILE_DEFINITIONS C4CORE_NO_FAST_FLOAT
    )
    add_dependencies(ryml_import ryml_build)

    add_library(ryml INTERFACE)
    add_dependencies(ryml ryml_build)
    target_include_directories(ryml INTERFACE ${RYML_INSTALL_DIR}/include)
    target_link_libraries(ryml INTERFACE ${RYML_INSTALL_DIR}/lib/libryml.a)
    target_compile_definitions(ryml INTERFACE C4CORE_NO_FAST_FLOAT)

    ExternalProject_Add_Step(ryml_build cleanup
      COMMAND ${CMAKE_COMMAND} -E rm -rf <SOURCE_DIR> <BINARY_DIR> <STAMP_DIR> <TMP_DIR> <DOWNLOAD_DIR>
      COMMENT "Cleaning ryml external project work directories"
      DEPENDEES install
      ALWAYS TRUE
    )
  else()
    # Everything is compiled on host; fetch rapidyaml directly
    message(STATUS "Fetching rapidyaml...")
    include(FetchContent)
    
    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON CACHE BOOL "" FORCE) # Suppress developer warnings from ryml
    
    FetchContent_Declare(rapidyaml
      GIT_REPOSITORY https://github.com/biojppm/rapidyaml.git
      GIT_TAG 47ec2fa184209687c20fd5bc05621e1cb1200311
    )
    FetchContent_MakeAvailable(rapidyaml)

    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS OFF CACHE BOOL "" FORCE)
  endif()

else()
  if(TARGET ryml::ryml)
    message(STATUS "Using ryml from upstream project.")
  else()
    message(FATAL_ERROR "RYML_IN_MIRCO=OFF but no ryml target was available in the upstream project.")
  endif()

endif()

add_library(mirco_utils
  src/mirco_utils.cpp
  )

target_link_libraries(mirco_utils PUBLIC ryml)

add_library(mirco_inputparameters_yaml
  src/mirco_inputparameters_yaml.cpp
  )
target_link_libraries(mirco_inputparameters_yaml PRIVATE mirco_inputparameters mirco_utils Kokkos::kokkos)

add_library(mirco_lib INTERFACE)
target_link_libraries(mirco_lib INTERFACE mirco_needKK mirco_core mirco_topology mirco_inputparameters)

target_include_directories(mirco_lib INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/>
  $<INSTALL_INTERFACE:include/mirco/>
  )

add_library(mirco::mirco_lib ALIAS mirco_lib)

target_compile_definitions(mirco_core PUBLIC
  $<$<BOOL:${MIRCO_ENABLE_VISUALIZATIONEXPORT}>:MIRCO_ENABLE_VISUALIZATIONEXPORT=1>
  $<$<NOT:$<BOOL:${MIRCO_ENABLE_VISUALIZATIONEXPORT}>>:MIRCO_ENABLE_VISUALIZATIONEXPORT=0>
)

# Compile mirco itself
add_executable(mirco src/main.cpp)
target_link_libraries(mirco PUBLIC mirco::mirco_lib mirco_inputparameters_yaml mirco_utils Kokkos::kokkos)

# Compile the flatMirco utility
if(MIRCO_ENABLE_FLAT)
  add_executable(flatMirco src/main_flat.cpp)
  target_link_libraries(flatMirco PUBLIC mirco_core mirco_topology mirco_shapefactors Kokkos::kokkos KokkosKernels::kokkoskernels)
endif()

# Install mirco (to be used as a library by other codes)
install(TARGETS mirco_lib mirco_needKK mirco_core mirco_topology mirco_inputparameters mirco_shapefactors
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
  target_link_libraries(tests PUBLIC mirco::mirco_lib mirco_inputparameters_yaml mirco_utils mirco_shapefactors gtest Kokkos::kokkos)

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
