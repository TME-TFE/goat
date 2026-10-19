# Configuration file for the grid deformation makefile. Variables etc 
# are listed below. When introducing new variables or other 
# functionality, please document using the '##' prefix. THis should
# be compatible with the 'help' target in the makefile. 

##
## %===================================================================%
## %                                                                   %
## %                             VARIABLES                             %
## %                                                                   %
## %===================================================================%
##
## Variables are defined in config.mk

# Target-specific variables
## COMPDIRVARS:        define variables for compiler directives
COMPDIRVARS = 
ifdef SOLPSTOP
COMPDIRVARS += -DSOLPS
endif

## % Library paths
## %==============
## LAPACKPATH 			: LAPACK library path (user defined)
LAPACKPATH ?= -lopenblas

## BLASPATH 			: BLAS library path (user defined)
BLASPATH ?= -lopenblas

## UMFPACKPATH 			: UMFPACK library path (user defined)
UMFPACKPATH ?= -lumfpack

## CXSparsePATH          : CXSparse library path (user defined)
CXSparsePATH ?=

## DMUMPSPATH            : DMUMPS library path (user defined, optional)
#DMUMPSPATH = -ldmumps -lmumps_common -L/usr/lib -lmetis  -lesmumps -L../../PORD/lib/ -lpord -L../../lib -L../../libseq -lscotch -lscotcherr -L../../libseq/libmpiseq.a -lpthread
#DMUMPSPATH = -ldmumps -lmumps_common -L/usr/lib -lmetis  -lesmumps -lpord -lscotch -lscotcherr -lpthread
#DMUMPSPATH = -L../../MUMPS_5.8.0/lib -lsmumps -ldmumps -lmumps_common -L/usr/lib  -lparmetis -lmetis -L../../PORD/lib/ -lpord -L/usr/lib -lptesmumps -lptscotch -lptscotcherr -lscalapack-openmpi  -llapack   -lblas -lpthread


## DMUMPSLIBPATH        : DMUMPS include path (user defined, optional)
ifdef DMUMPS_LPATH 
ifdef DMUMPS_IPATH
ifndef NO_USE_MPI
            # Define mumps 
            MUMPS = yes 
            # Mumps has to be compiled with MPI
            USE_MPI = yes 
            $(info % MUMPS library and include paths set, MPI available. Compiling with MUMPS) 
else 
USE_MPI =
DMUMPS_LPATH =
DMUMPS_IPATH =
            $(info % MUMPS paths available, but no MPI. Not compiling MUMPS.)
endif
else
        # Not all paths define, issue message and undefine to ensure proper compilation
        $(info % MUMPS include path not set, set "DMUMPS_LPATH" and "DMUMPS_IPATH" to enable compilation with MUMPS)  
DMUMPS_LPATH =
DMUMPS_IPATH =
endif
else
    # Not all paths define, issue message and undefine to ensure proper compilation
    $(info % MUMPS library path not set, set "DMUMPS_LPATH" and "DMUMPS_IPATH" to enable compilation with MUMPS) 
DMUMPS_LPATH =
DMUMPS_IPATH =
endif

# Define MUMPS for the compiler
ifdef MUMPS 
COMPDIRVARS += -DMUMPS 
endif

## SOLPSTOP            : path to SOLPS (overridden if SOLPSTOP is define)
ifdef SOLPSTOP
DOSOLPS = true
endif
ifdef HOST_NAME
else
HOST_NAME = DEFAULT
endif 

# Default prefix for OBJDIR: standalone
PREF_OBJDIR = standalone
ifdef USE_EIRENE
PREF_OBJDIR = couple_SOLPS-ITER
endif

# Extensions for SOLPS object directories when various options are used
ifdef USE_IMPGYRO
EXT_IMPGYRO = .ig
else
ifdef USE_MPI
EXT_MPI = .mpi
endif
endif
ifdef USE_OPENMP
EXT_OPENMP = .openmp
endif
ifdef SOLPS_DEBUG
GOAT_DEBUG = yes
EXT_DEBUG = .debug
IMAS_AMNS_DEBUG = yes
else
IMAS_AMNS_DEBUG = no
endif
ifdef DIFF_D
EXT_DIFF = .diff_d
DIFF = yes
DIFFDIR = builds/differentiated_files${EXT_DIFF}
endif
ifdef DIFF_B
EXT_DIFF = .diff_b
DIFF = yes
DIFFDIR = builds/differentiated_files${EXT_DIFF}
endif
ifdef TGT
EXT_DIFF = .tgt
DIFF = yes
DIFFDIR = src/differentiation/tangent
endif
ifdef ADJ
EXT_DIFF = .adj
DIFF = yes
DIFFDIR = src/differentiation/adjoint
endif
ifdef ADJ_SHAPE
EXT_DIFF = .adj_shape
DIFF = yes
DIFFDIR = src/differentiation/adjoint_shape
endif
ifeq ($(strip $(GOAT_DEBUG)),yes)
EXT_DEBUG = .debug
endif
##
## % Compiler
## %=========
## FC			: Compiler to be used for fortran (overridden if COMPILER is defined)
ifdef COMPILER
ifeq ($(strip $(COMPILER)),gfortran)
ifdef USE_MPI
            FC = mpif90 
else
            FC = gfortran 
endif  
else ifeq ($(strip $(COMPILER)),ifort64)
ifdef USE_MPI
            FC = mpiifort
else 
            FC = ifort 
endif 
else
ifdef USE_MPI
            FC = mpif90 
else
            FC = gfortran 
endif  
endif 
else 
    COMPILER = gfortran
ifdef USE_MPI
        FC = mpif90 
else
        FC = gfortran 
endif  
endif 

# Directory where objectcode/binaries will be created
B25LIBPATH = ${SOLPSTOP}/modules/B2.5/builds/${PREF_OBJDIR}.${HOST_NAME}.${COMPILER}${EXT_OPENMP}${EXT_MPI}${EXT_IMPGYRO}${EXT_DIFF}${EXT_DEBUG}

# Build path
BUILDDIR = ${PREF_OBJDIR}.${HOST_NAME}.${COMPILER}${EXT_OPENMP}${EXT_MPI}${EXT_IMPGYRO}${EXT_DIFF}${EXT_DEBUG}

## % Include paths
## %==============
## SUITESPARSEPATH      : SuiteSparse header file path
SUITESPARSEPATH ?= -I/usr/include/suitesparse

## CFLAGS			: Compiler flags for standard compilation (may be overridden)
CFLAGS_DEF = -c -fopenmp
CFLAGS_DEF_NO_OMP = -c -g -Wall -O0 -Wno-unused-dummy-argument -Wno-maybe-uninitialized -fcheck=all -Wno-uninitialized 
## CFLAGS_OMP	: compiler flags for OpenMP 
CFLAGS_OMP = -c -Wall -fopenmp
## CFLAGS_DEBUG		: compiler flags for debugging
CFLAGS_DEBUG = -c -g -Wall -O0 -Wno-unused-dummy-argument -Wno-maybe-uninitialized -fcheck=all -fopenmp -Wno-uninitialized
## CFLAGS_OMP_DEBUG	: compiler flags for OpenMP and debugging
CFLAGS_OMP_DEBUG = -c -g -Wall  -O0 -fopenmp
CFLAGS_PERF = -c -O2 -Wno-unused-dummy-argument -Wno-maybe-uninitialized -fopenmp -Wno-uninitialized

## CC           : Compiler to be used for C
ifdef USE_MPI
    CC = mpicc
else 
    CC = gcc 
endif 

CCFLAGS_DEF = -c -g
CCFLAGS_PERF = -c -Wall -O2

## % Linker
## %=======
## LFLAGS			: linking flags to be used (apart from libraries)
LFLAGS_DEF =   -fcheck=all -fopenmp
LFLAGS_DEF_NO_OMP =   -fcheck=all
## LFLAGS_DEBUG 	: linking flags for debugging
LFLAGS_DEBUG = -g -fcheck=all

## LFLAGS			: linking flags to be used for openMP
LFLAGS_OMP = -g -fopenmp
LFLAGS_PERF = -fopenmp

## LFLAGS_DEBUG 	: linking flags for debugging
LFLAGS_OMP_DEBUG = -g -fopenmp

# Set flags
#==========
# Set the CFLAGS
CFLAGS = $(CFLAGS_DEF) $(COMPDIRVARS)
ifeq ($(strip $(GOAT_DEBUG)),yes)
CFLAGS += -g -O0 -Wall -Wno-unused-dummy-argument -Wno-maybe-uninitialized -fcheck=all -Wno-uninitialized
else
CFLAGS += -O3
endif

# Set the linking flags
LFLAGS = $(LFLAGS_DEF)

# Set CFLAGS for C compiler
CCFLAGS = $(CCFLAGS_DEF) $(COMPDIRVARS)
ifeq ($(strip $(GOAT_DEBUG)),yes)
CCFLAGS += -g -Wall -O0 
else
CCFLAGS += -O3
endif

##
## % Files
## %======
## GENERAL_FILES				: All general files (e.g. precision definition, ... )
GENERAL_FILES = src/General/mod_errorhandler.F90 src/General/mod_readwrite.F90 src/General/mod_inputfileparser.F90 src/General/mod_plotter.F90 \
    src/General/mod_sparseinterface.F90  src/General/mod_sort.F90 $(sort $(wildcard src/General/*.F90))      

## DRIVER_FILES			: Driver filenames (.F90) - unsequenced
DRIVER_FILES = $(sort $(wildcard src/Drivers/Goat/*.F90))

## SODRIVER_FILES			: Shape optimization driver filenames (.F90) - unsequenced
SODRIVER_FILES = $(sort $(wildcard src/Drivers/ShapeOpt/*.F90))

## MODULE_FILES			: Module filenames (.F90, .F) - sequence matters
MODULE_FILES_GOAT =  src/Modules/Goat/goatmod_userinput.F90 src/Modules/Goat/goatmod_types.F90\
    $(sort $(wildcard src/Modules/Goat/*.F90))
MODULE_FILES_GD = src/Modules/GD/gdmod_types.F90 src/Modules/GD/gdmod_userinput.F90 src/Modules/GD/gdmod_plots.F90 src/Modules/GD/gdmod_designvariables.F90 \
    src/Modules/GD/gdmod_utility_optimization.F90 src/Modules/GD/gdmod_constraints.F90\
    $(sort $(wildcard src/Modules/GD/*.F90))
MODULE_FILES_GA = src/Modules/GA/gamod_math.F90 src/Modules/GA/gamod_types.F90 src/Modules/GA/gamod_utility.F90 src/Modules/GA/gamod_driver.F90  \
    $(sort $(wildcard src/Modules/GA/*.F90))    
MODULE_FILES_GG = src/Modules/GG/ggmod_topology2D.F90 src/Modules/GG/ggmod_vertexdistribution2D.F90 \
    src/Modules/GG/ggmod_solpsregions.F90 \
    src/Modules/GG/ggmod_gridgeneration2D.F90 src/Modules/GG/ggmod_gridgenerator.F90

MODULE_FILES = $(sort $(wildcard src/Modules/Goat/*.F90))\
    src/Modules/GD/gdmod_types.F90 src/Modules/GD/gdmod_userinput.F90 src/Modules/GD/gdmod_plots.F90 src/Modules/GD/gdmod_designvariables.F90 \
    src/Modules/GD/gdmod_utility_optimization.F90 src/Modules/GD/gdmod_constraints.F90\
    $(sort $(wildcard src/Modules/GD/*.F90)) \
    $(sort $(wildcard src/Modules/GA/*.F90)) \
    $(sort $(wildcard src/Modules/*.F90)) $(sort $(wildcard src/Modules/*.F)) \
    src/Modules/GG/ggmod_topology2D.F90 src/Modules/GG/ggmod_vertexdistribution2D.F90 \
    src/Modules/GG/ggmod_solpsregions.F90 \
    src/Modules/GG/ggmod_gridgeneration2D.F90

## AUXILIARY_FILES			: Auxiliary filenames (.F90) - unsequenced
AUXILIARY_FILES =  src/Auxiliary/mod_structured2Dgridding.F90 \
    $(sort $(wildcard src/Auxiliary/*.F90)) \
    src/Auxiliary/Interpolation/Interpolant2D_auxiliaries.F90 \
    src/Auxiliary/Interpolation/Interpolant2D.F90 \
    src/Auxiliary/Interpolation/Interpolant1D.F90 \
    src/Auxiliary/Interpolation/BicubicSplineInterpolant.F90 \
    src/Auxiliary/Interpolation/StructuredInterpolant2D.F90 \
    $(sort $(wildcard src/Auxiliary/Interpolation/*.F90)) \
    $(sort $(wildcard src/Auxiliary/Contour/*.F90)) \
    src/Auxiliary/mod_streamlinetracing2D.F90 \
    src/Auxiliary/Graphs/mod_graph.F90

## SETUP_FILES			: setup file generation names (.F90) - unsequenced
SETUP_FILES = $(sort $(wildcard src/Setup/*.F90))

## NUMERICS_FILES 		: numeric file generation names (.F90) - unsequenced
NUMERICS_FILES = src/Numerics/PolygonLevelsetFunction2D.F90 $(sort $(wildcard src/Numerics/*.F90))

## OPTIMIZATION 		: optimization file generation names (.F90) - unsequenced
OPTIMIZATION_FILES = src/Optimization/optmod_designvariables.F90 src/Optimization/optmod_costfunction.F90 \
    src/Optimization/optmod_constraints.F90 src/Optimization/optmod_monitor.F90 \
    src/Optimization/optmod_numerics.F90 $(sort $(wildcard src/Optimization/*.F90))

## CONSTANTS            : constants such as precision and special characters (.F90) - unsequenced
CONSTANTS_FILES = src/Constants/mod_global_environment.F90 src/Constants/mod_precision.F90 $(sort $(wildcard src/Constants/*.F90))

## Clayer               : c files for interfacing with other c code
CLAYER_FILES    = $(sort $(wildcard src/Clayer/*.c))

## ClayerF              : fortran files for interfacing with other c code
CLAYERF_FILES    =   src/Clayer/CUtilities.F90 src/Clayer/CSparseF.F90 src/Clayer/Clayer.F90

## ShapeOpt             : fortran files for shape optimization
SHAPEOPT_FILES  = src/Modules/ShapeOpt/somod_userinput.F90 \
    src/Modules/ShapeOpt/somod_designvariables.F90 src/Modules/ShapeOpt/somod_costfunction.F90 \
    src/Modules/ShapeOpt/somod_constraints.F90   src/Modules/ShapeOpt/somod_optimizationengine.F90

## ShapeOptSolps            : fortran files for shape optimization with SOLPS
SHAPEOPTSOLPS_FILES  =  src/Modules/ShapeOpt/somod_userinput.F90 \
    src/Modules/ShapeOpt/somod_designvariables.F90 src/Modules/ShapeOpt/somod_costfunction.F90 \
    src/Modules/ShapeOpt/somod_constraints.F90  src/Modules/ShapeOpt/sosmod_costfunction.F90 \
    src/Modules/ShapeOpt/somod_optimizationengine.F90 
    

## % Targets
## %========
## GOAT_TARGETS             : Targets to be run for the full goat
GOAT_TARGETS = Clayer ClayerF Constants General Auxiliary Numerics Optimization Modules  \
    Drivers 

## TEST_TARGETS             : Targets to be run for goat tests
TEST_TARGETS = $(GOAT_TARGETS) 

## CTEST_TARGETS            : Targets to be run to test C layer
CTEST_TARGETS = Clayer

## SHAPEOPT_TARGETS         : Targets to be run for shape optimization program
ifdef DOSOLPS
SHAPEOPT_TARGETS = Clayer ClayerF Constants General Auxiliary Numerics Optimization Modules  \
    ShapeOptimizationSolps Drivers SODrivers
else
SHAPEOPT_TARGETS = Clayer ClayerF Constants General Auxiliary Numerics Optimization Modules  \
    ShapeOptimization Drivers SODrivers
endif
