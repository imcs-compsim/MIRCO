# parameters
set(CTEST_SOURCE_DIRECTORY "/home_local/github-runner/actions-runner/_work/MIRCO/MIRCO")
set(CTEST_BINARY_DIRECTORY "/home_local/github-runner/actions-runner/_work/MIRCO/mirco_build")
set($ENV{LC_MESSAGES} "en_EN" ) # set output to english such that ctest can analyze it

set(CTEST_SITE "$ENV{HOSTNAME}")
set(CTEST_BUILD_NAME "$ENV{CTEST_BUILD_NAME_GITLAB}")

# prepare environment
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(WITH_MEMCHECK FALSE)
set(WITH_COVERAGE FALSE)

# ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

# prepare commands
set(CTEST_CONFIGURE_COMMAND "${CTEST_SOURCE_DIRECTORY}/do-configure")
set(CTEST_BUILD_COMMAND "$ENV{CTEST_MAKE}")
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS "100")

# do the testing
ctest_start("${HOSTNAME}")
ctest_configure()
ctest_build(RETURN_VALUE testBuild NUMBER_WARNINGS numWarnings)
ctest_test(RETURN_VALUE testRes)

if (NOT ${testBuild} EQUAL 0) # Send error for a failed build
  message( SEND_ERROR "Mirco build failed!" )
endif (NOT ${testBuild} EQUAL 0)

if (NOT ${numWarnings} EQUAL 0 AND $ENV{CTEST_FAIL_ON_WARNING} EQUAL 1) # Send error for warnings if enabled
  message( SEND_ERROR "Mirco build issued build warnings!" )
endif (NOT ${numWarnings} EQUAL 0 AND $ENV{CTEST_FAIL_ON_WARNING} EQUAL 1)

if (NOT ${testRes} EQUAL 0)
  message( SEND_ERROR "MIRCO tests failed!" )
endif (NOT ${testRes} EQUAL 0)

