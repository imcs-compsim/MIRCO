##############################################################################
#
# Make variables for use by Trilinos/ML clients. 
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################


## ---------------------------------------------------------------------------
## Compilers used by Trilinos/ML build
## ---------------------------------------------------------------------------

ML_CXX_COMPILER=/usr/bin/mpic++

ML_C_COMPILER=/usr/bin/mpicc

ML_FORTRAN_COMPILER=/usr/bin/mpif90
ML_Fortran_COMPILER=/usr/bin/mpif90


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/ML build
## ---------------------------------------------------------------------------

## Set compiler flags, including those determined by build type
ML_CXX_FLAGS= -pedantic -Wall -Wno-long-long -Wwrite-strings  -std=c++11 -O3
ML_CXX_COMPILER_FLAGS= -pedantic -Wall -Wno-long-long -Wwrite-strings  -std=c++11 -O3

ML_C_FLAGS= -pedantic -Wall -Wno-long-long -std=c99 -O3
ML_C_COMPILER_FLAGS= -pedantic -Wall -Wno-long-long -std=c99 -O3

ML_FORTRAN_FLAGS= -O3
ML_Fortran_COMPILER_FLAGS= -O3

## Extra link flags (e.g., specification of fortran libraries)
ML_EXTRA_LD_FLAGS=

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty. 
ML_SHARED_LIB_RPATH_COMMAND=
ML_BUILD_SHARED_LIBS=OFF

ML_LINKER=/usr/bin/ld
ML_AR=/usr/bin/ar

## ---------------------------------------------------------------------------
## Set library specifications and paths 
## ---------------------------------------------------------------------------

## The project version number
ML_VERSION=12.3

## List of package include dirs
ML_INCLUDE_DIRS= -I/rzhome/nas/compsim/public/lib/Q1_2015/TPL/trilinos-build/include

## List of package library paths
ML_LIBRARY_DIRS= -L/rzhome/nas/compsim/public/lib/Q1_2015/TPL/trilinos-build/lib

## List of package libraries
ML_LIBRARIES= -lml -lifpack -lamesos -laztecoo -lisorropia -lepetraext -ltpetraext -ltpetrainout -ltpetra -lkokkostsqr -ltpetrakernels -ltpetraclassiclinalg -ltpetraclassicnodeapi -ltpetraclassic -ltriutils -lzoltan -lepetra -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lkokkosalgorithms -lkokkoscontainers -lkokkoscore

## Specification of directories for TPL headers
ML_TPL_INCLUDE_DIRS= -I/rzhome/nas/compsim/public/lib/Q1_2015/TPL/include -I/rzhome/nas/compsim/public/lib/Q1_2015/TPL/include -I/rzhome/nas/compsim/public/lib/Q1_2015/TPL/include -I/usr/include

## Specification of directories for TPL libraries
ML_TPL_LIBRARY_DIRS=

## List of required TPLs
ML_TPL_LIBRARIES= /rzhome/nas/compsim/public/lib/Q1_2015/TPL/lib/libumfpack.a /rzhome/nas/compsim/public/lib/Q1_2015/TPL/lib/libamd.a /rzhome/nas/compsim/public/lib/Q1_2015/TPL/lib/libsuperlu_dist_2.5.a /rzhome/nas/compsim/public/lib/Q1_2015/TPL/lib/libparmetis.a /rzhome/nas/compsim/public/lib/Q1_2015/TPL/lib/libmetis.a /usr/lib/liblapack.so /usr/lib/libblas.so

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

ML_MPI_LIBRARIES=
ML_MPI_LIBRARY_DIRS=
ML_MPI_INCLUDE_DIRS=
ML_MPI_EXEC=/usr/bin/mpiexec
ML_MPI_EXEC_MAX_NUMPROCS=4
ML_MPI_EXEC_NUMPROCS_FLAG=-np

## ---------------------------------------------------------------------------
## Set useful general variables 
## ---------------------------------------------------------------------------

## The packages enabled for this project
ML_PACKAGE_LIST= ML Ifpack Amesos AztecOO Isorropia EpetraExt Tpetra TpetraCore TpetraTSQR TpetraKernels TpetraClassic Triutils Zoltan Epetra Teuchos TeuchosKokkosComm TeuchosKokkosCompat TeuchosRemainder TeuchosNumerics TeuchosComm TeuchosParameterList TeuchosCore KokkosAlgorithms KokkosContainers KokkosCore

## The TPLs enabled for this project
ML_TPL_LIST= UMFPACK SuperLUDist ParMETIS Boost LAPACK BLAS MPI Pthread

## ---------------------------------------------------------------------------
## Deprecated variables
##   These variables from the autotools version of Makefile.export.<package>
##   have been deprecated. To help avoid confusion with projects that haven't
##   migrated to the new version of the variable names when the don't compile.
##   This will put an obvious statement in the link line that these Trilinos
##   Makefile variables are no longer used. Without this there would be a link
##   error with no obvious reason as to why other than the lack of Trilinos libs
##   which would be confusing since the link line in the Makefile would have 
##   stated to put them on the line.
## ---------------------------------------------------------------------------

ML_LIBS="The Trilinos variable ML_LIBS has been deprecated. Please see the documentation on the Makefile.export capability for instructions on how to transition to the new system."
ML_INCLUDES="The Trilinos variable ML_INCLUDES has been deprecated. Please see the documentation on the Makefile.export capability for instructions on how to transition to the new system."


