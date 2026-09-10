# A very simple and straightforward makefile without too many options 
# yet. To be improved in the future. An attempt has been made to clean
# up the compilation and linking steps etc by moving all .o and .mod
# files into a build directory that is constructed during compilation

# Preamble
#=========
# Set the executable name
EXEC_NAME = $(MAKECMDGOALS).exe

# Check if we are debugging
$(info % ===========================)
$(info % Makefile preamble execution)
$(info % ===========================)
$(info % Makefile command: $(MAKECMDGOALS)) 

# Check debug mode
ifeq ($(findstring _debug, $(MAKECMDGOALS)), _debug)
GOAT_DEBUG = yes 
else
GOAT_DEBUG = no 
endif
$(info % Goat debug mode: $(GOAT_DEBUG))

# Check if MUMPS is set
ifdef MUMPS 
$(info % Compiling with MUMPS solver)
endif

# Include the config file
include config.mk

# Echo the different environment variables 
$(info % ===========================)
$(info % Executing the goat makefile)
$(info % ===========================)
$(info % )
$(info % General environment variables)
$(info % -----------------------------)
$(info % fortran compiler: $(FC))
$(info % C compiler: $(CC))
ifdef USE_MPI 
$(info % Running with MPI)
endif
$(info % Build directory: $(BUILDDIR))
$(info % )
$(info % SOLPS specific environment variables)
$(info % ------------------------------------)
$(info % running at host: $(HOST_NAME))
$(info % solps path: $(SOLPSTOP))
$(info % B2.5 build path: $(BUILDDIR))
$(info % B2.5 libary path: $(B25LIBPATH))

## %===================================================================%
## %                                                                   %
## %                          PRELIMINARIES                            %
## %                                                                   %
## %===================================================================%
## Construction of directories
## %==========================
## CREATE_BUILDDIR: construct build directory based on BUILDDIR environment variable
CREATE_BUILDDIR:= $(shell mkdir ./builds)
CREATE_BUILDDIR:= $(shell mkdir ./builds/$(BUILDDIR))

## CREATE_EXEDIR: construct executable directory 
CREATE_EXEDIR:= $(shell mkdir ./executables)

BUILDDIR :=./builds/$(BUILDDIR)


##
## %===================================================================%
## %                                                                   %
## %                             TARGETS                               %
## %                                                                   %
## %===================================================================%
##
## % Main targets
## %=============
# Linking
## goat			: Create executable for goat
## tests		: Create executable for tests
## testc 		: Create executable for tests of C layer
## goattranslator: Create executable for GOAToptions file translator
## shapeopt 	: Create executable for shape optimization with goat

goat: $(addprefix $(BUILDDIR)/, $(GOAT_TARGETS)) $(BUILDDIR)/goat.o
	-mv -f *.o $(BUILDDIR);  
	-mv -f *.mod $(BUILDDIR); 
	$(FC) $(LFLAGS) -o $(BUILDDIR)/$(EXEC_NAME) $(BUILDDIR)/*.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) $(DMUMPS_LPATH) $(CXSparsePATH) -lcxsparse \
	$(SUITESPARSEPATH) -I src/Clayer/Include $(DMUMPS_IPATH); 
	rm $(BUILDDIR)/Goat.o; 
	cp $(BUILDDIR)/$(EXEC_NAME) ./executables/.

goat_debug: goat

tests: $(addprefix $(BUILDDIR)/,$(TEST_TARGETS) ) $(BUILDDIR)/tests.o 
	-mv -f *.o $(BUILDDIR);  
	-mv -f *.mod $(BUILDDIR); 
	$(FC) $(LFLAGS) -o $(BUILDDIR)/tests.exe $(BUILDDIR)/*.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) $(DMUMPS_LPATH) $(CXSparsePATH) -lcxsparse \
	$(SUITESPARSEPATH) -I src/Clayer/Include
	rm $(BUILDDIR)/tests.o; 
	cp $(BUILDDIR)/tests.exe ./executables/.

testc: $(addprefix $(BUILDDIR)/,$(CTEST_TARGETS) ) $(BUILDDIR)/testc.o 
	-mv -f *.o $(BUILDDIR);  
	-mv -f *.mod $(BUILDDIR); 
	$(CC) -o $(BUILDDIR)/testc.exe $(BUILDDIR)/*.o $(CXSparsePATH) -lcxsparse $(SUITESPARSEPATH) -I src/Clayer/Include
	rm $(BUILDDIR)/Testc.o; 
	cp $(BUILDDIR)/testc.exe ./executables/.

shapeopt: $(addprefix $(BUILDDIR)/, $(SHAPEOPT_TARGETS) ) $(BUILDDIR)/shapeopt.o 
	-mv -f *.o $(BUILDDIR);  
	-mv -f *.mod $(BUILDDIR); 
ifdef DOSOLPS
	$(FC) $(LFLAGS) -o $(BUILDDIR)/$(EXEC_NAME) $(BUILDDIR)/*.o $(B25LIBPATH)/adStack.o \
	 $(B25LIBPATH)/b2mod_cdf.o $(B25LIBPATH)/smax.o $(B25LIBPATH)/smin.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) $(DMUMPS_LPATH) \
	 $(CXSparsePATH) -lcxsparse $(SUITESPARSEPATH) -I src/Clayer/Include  -I$(B25LIBPATH) -L$(B25LIBPATH) -l:libb2.a -L$(B25LIBPATH) -l:libb2.a -lnetcdf $(LD_NETCDF)
else
	$(FC) $(LFLAGS) -o $(BUILDDIR)/$(EXEC_NAME) $(BUILDDIR)/*.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) $(DMUMPS_LPATH) $(CXSparsePATH) -lcxsparse \
	$(SUITESPARSEPATH) -I src/Clayer/Include
endif
	rm $(BUILDDIR)/ShapeOptimization.o; 
	cp $(BUILDDIR)/$(EXEC_NAME) ./executables/.

shapeopt_debug: shapeopt

#shapeopt_solps: $(addprefix $(BUILDDIR)/,$(SHAPEOPTSOLPS_TARGETS) ) $(BUILDDIR)/shapeopt_solps.o 
#	-mv -f *.o $(BUILDDIR);  
#	-mv -f *.mod $(BUILDDIR); 
#	$(FC) $(LFLAGS) -o $(BUILDDIR)/shapeopt_solps $(BUILDDIR)/*.o $(B25LIBBPATH)/adStack.o \
	 $(B25LIBBPATH)/b2mod_cdf.o $(LAPACKPATH) $(BLASPATH) $(UMFPACKPATH) \
	 $(CXSparsePATH) -lcxsparse $(SUITESPARSEPATH) -I src/Clayer/Include  -I$(B25LIBPATH) -L$(B25LIBPATH) -l:libb2.a -L$(B25LIBPATH) -l:libb2.a -lnetcdf $(LD_NETCDF)
#	rm $(BUILDDIR)/ShapeOptimization.o; 
#	cp $(BUILDDIR)/shapeopt_solps.exe ./executables/.


## % Runfiles
## %=========
## Goat.o 			: main run file goat
$(BUILDDIR)/goat.o: Runfiles/Goat.F90
	$(FC) $(CFLAGS) Runfiles/Goat.F90 -I$(BUILDDIR) $(DMUMPS_IPATH) $(DMUMPS_LPATH)

## Tests.o 			: all tests
$(BUILDDIR)/tests.o: Runfiles/Tests.F90 
	$(FC) $(CFLAGS) Runfiles/Tests.F90 -I$(BUILDDIR) $(DMUMPS_IPATH)

## Testc.o 			: C layer tests
$(BUILDDIR)/testc.o: Runfiles/Testc.c 
	$(CC) $(CCFLAGS) Runfiles/Testc.c $(SUITESPARSEPATH) -I src/Clayer/Include -I$(BUILDDIR) $(DMUMPS_IPATH)

## shapeopt.o		: shape optimization 
$(BUILDDIR)/shapeopt.o : Runfiles/ShapeOptimization.F90 
	$(FC) $(CFLAGS) Runfiles/ShapeOptimization.F90 -I$(BUILDDIR) $(DMUMPS_IPATH)

## shapeopt_solps.o		: shape optimization with solps
$(BUILDDIR)/shapeopt_solps.o : Runfiles/ShapeOptimization.F90 
	$(FC) $(CFLAGS) Runfiles/ShapeOptimization.F90 -I$(BUILDDIR) $(DMUMPS_IPATH)

##
## % Folder compilation targets
## %===========================
## Constants 		: compile files containing constants
$(BUILDDIR)/Constants: $(CONSTANTS_FILES) 
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR) 
	touch $(BUILDDIR)/Constants 

## General 			: compile general files and modules
$(BUILDDIR)/General: $(GENERAL_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR) $(DMUMPS_IPATH)
	touch $(BUILDDIR)/General
	
## Modules			: compile modules
$(BUILDDIR)/Modules_goat: $(MODULE_FILES_GOAT)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR) 
	touch $(BUILDDIR)/Modules_goat
$(BUILDDIR)/Modules_GD: $(MODULE_FILES_GD)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR) 
	touch $(BUILDDIR)/Modules_GD
$(BUILDDIR)/Modules_GA: $(MODULE_FILES_GA)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR) 
	touch $(BUILDDIR)/Modules_GA	
$(BUILDDIR)/Modules_GG: $(MODULE_FILES_GG)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR) 
	touch $(BUILDDIR)/Modules_GG
$(BUILDDIR)/Modules: $(BUILDDIR)/Modules_goat $(BUILDDIR)/Modules_GD $(BUILDDIR)/Modules_GA\
	$(BUILDDIR)/Modules_GG 

## Auxiliary			: compile auxiliary routines
$(BUILDDIR)/Auxiliary: $(AUXILIARY_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR) 
	touch $(BUILDDIR)/Auxiliary

## Drivers			: compile all (goat) driver routines
$(BUILDDIR)/Drivers: $(DRIVER_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR) 
	touch $(BUILDDIR)/Drivers

## SODrivers			: compile all shape optimization driver routines
$(BUILDDIR)/SODrivers: $(SODRIVER_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR) 
	touch $(BUILDDIR)/SODrivers

## Setup			: compile all setup routines
$(BUILDDIR)/Setup: $(SETUP_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR) 
	touch $(BUILDDIR)/Setup

## Numerics		 	: compile numerics routines
$(BUILDDIR)/Numerics: $(NUMERICS_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)  
	touch $(BUILDDIR)/Numerics

## Optimization		 	: compile optimization routines
$(BUILDDIR)/Optimization: $(OPTIMIZATION_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR) 
	touch $(BUILDDIR)/Optimization

## Clayer 				: compile C interlayer routines
$(BUILDDIR)/Clayer: $(CLAYER_FILES)
	$(CC) $(CCFLAGS) $^ $(SUITESPARSEPATH) -I src/Clayer/Include -I$(BUILDDIR)
	touch $(BUILDDIR)/Clayer

## ClayerF 				: compile fortran modules in Clayer
$(BUILDDIR)/ClayerF: $(CLAYERF_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/ClayerF

## ShapeOptimization 			: compile shape optimization modules
$(BUILDDIR)/ShapeOptimization: $(SHAPEOPT_FILES)
	$(FC) $(CFLAGS) $^ -I$(BUILDDIR)
	touch $(BUILDDIR)/ShapeOptimization

## ShapeOptimizationSolps 			: compile shape optimization modules for SOLPS
$(BUILDDIR)/ShapeOptimizationSolps: $(SHAPEOPTSOLPS_FILES)
	$(FC) $(CFLAGS) -I$(B25LIBPATH) -I$(BUILDDIR) -L$(B25LIBPATH) -l:libb2.a -lnetcdf $^ 
	touch $(BUILDDIR)/ShapeOptimizationSolps 



##
## % Run commands
## %=============
## gd			: Run grid deformation with gdrun
# Run
.PHONY: gd
gd: 
	gdrun
	make clean

##
## % Auxiliary targets
## %==================
## clean			: clean by removing *.o 
# Cleanup
.PHONY: clean
clean: 
	rm -f *.o; 
	find . -name "*.mod" -type f -delete; \
	rm -f $(GOAT_TARGETS); \
	rm -f $(SHAPEOPT_TARGETS); \
	rm -f $(SHAPEOPTSOLPS_TARGETS); \
	rm -f goat; \
	rm -f tests; \
	rm -f testc; \
	rm -f shapeopt; \

## cleanbuilds   		: clean builds directory (removes all subdirectories)
.PHONY: cleanbuilds 
cleanbuilds: 
	rm -rf ./builds/*

## deepclean   		    : execute clean and remove all build directories
.PHONY: deepclean 
deepclean: cleanbuilds clean 

	
## help			: print out documentation
# Help - prints out all the ## statements
.PHONY : help 
help : config.mk Makefile  
	@sed -n 's/^##//p' $^

##
##
