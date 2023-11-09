print-%: ; @echo $*=$($*)
.PHONY: tests all clean
.DEFAULT_GOAL = tests

MODELERC ?= $(HOMEDIR)/.modelErc
ifneq ($(MPI),YES)
  USE_MPI=NO
else
  USE_MPI=YES
endif
FFLAGS = $(CPPFLAGS) $(F90FLAGS)

ifeq ($(wildcard $(MODELERC)),)
   $(error $(MODELERC) not found.)
endif
include $(MODELERC)
ifeq ($(COMPILER),)
   $(error you need to specify COMPILER in your ~.modelErc)
endif
SOURCE_DIR = ../model

# Unfortunately COMPILER is defined differently in modelE (gfortran)
# than in PFUNIT (GNU). So, we rename here to inherit the correct PFUNIT
# configuration.
ifeq ($(COMPILER),gfortran)
   COMPILER=GNU
endif

COMPLIBS = $(shell perl -e 'print reverse <>;' < $(SOURCE_DIR)/.liblist)
LIBS := $(foreach ITEM, $(COMPLIBS), $(SOURCE_DIR)/$(ITEM))
LIBS += $(SOURCE_DIR)/libmodel.a

# Determine modelE components to be tested based on availability
ifeq ($(COMPONENTS),)
  TESTDIRS = shared profiler MPI_Support
else
  TESTDIRS := $(foreach ITEM, $(COMPONENTS), $(wildcard $(ITEM)))
endif

TESTLIBS := $(foreach ITEM, $(TESTDIRS), $(ITEM)/lib$(ITEM).a)
TESTLIBS += model/libmodel.a
# The include paths are different when compiling tests and executable
TESTINCS := $(foreach ITEM, $(TESTDIRS), -I../$(SOURCE_DIR)/$(ITEM))
TESTINCS += -I../../model

FFLAGS += $(TESTINCS) -I$(PFUNIT)/mod -I$(SOURCE_DIR)
EXEFLAGS += $(FFLAGS) -I. -I../model/shared
EXEFLAGS += -DPFUNIT_EXTRA_INITIALIZE=initModelEtests

# Required libraries
PFLIB = -L$(PFUNIT)/lib -lpfunit
IOLIB += -L$(NETCDFHOME)/lib -lnetcdf 

ifeq ($(COMPILER),intel)
  FFLAGS += -DPROTECT_PARENS -O0 -g -traceback
else
  FFLAGS += -DPROTECT_PARENS -O0 -g -fbacktrace 
endif

# For MPI, use compiler wrapper and link to parallel-netcdf library
ifeq ($(USE_MPI),YES)
  FC = mpif90
  IOLIB += -L$(PNETCDFHOME)/lib -lpnetcdf
  EXEFLAGS += -DUSE_MPI
endif

# Targets:

EXE = tests.x

MYLIBS = $(TESTLIBS) $(LIBS) $(PFLIB) $(IOLIB)
$(EXE): testSuites.inc TESTS
	$(FC) -o $@ initModelEtests.F90 $(EXEFLAGS) \
	$(PFUNIT)/include/driver.F90 $(MYLIBS)

all: $(EXE)
ifeq ($(USE_MPI),YES)
	mpirun -np 12 ./$(EXE)
else
	./$(EXE)
endif

TESTS:
	@for dir in $(TESTDIRS); do \
	make -C $$dir all; done
	make -C model all

tests: all

ALLDIRS = profiler shared MPI_Support tracers model
clean:
	@for dir in $(ALLDIRS); do \
	make -C $$dir clean; done
	$(RM) tests.x

export FC
export MPIF90
export FFLAGS
export USE_MPI
