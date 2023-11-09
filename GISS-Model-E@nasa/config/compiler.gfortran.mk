
ifneq ($(MODELE_CPATH),)
ifneq ($(CPATH),)
  CPATH_HACK = CPATH=$(strip $(MODELE_CPATH)):$(CPATH)
else
  CPATH_HACK = CPATH=$(strip $(MODELE_CPATH))
endif
endif

ifneq ($(MODELE_LIBRARY_PATH),)
ifneq ($(LIBRARY_PATH),)
  LIBRARY_PATH_HACK = LIBRARY_PATH=$(strip $(MODELE_LIBRARY_PATH)):$(LIBRARY_PATH)
else
  LIBRARY_PATH_HACK = LIBRARY_PATH=$(strip $(MODELE_LIBRARY_PATH))
endif
endif

# LIBRARY_PATH is specified for F90 because it is used for linking
F90 = $(LIBRARY_PATH_HACK) gfortran
CC = $(CPATH_HACK) gcc
ifneq ($(CPATH_HACK),)
  CPP := $(CPATH_HACK) $(CPP)
endif

FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CPPFLAGS += -DCOMPILER_G95
FFLAGS = -g -cpp -fconvert=big-endian -O2 -fno-range-check
F90FLAGS = -g -cpp -fconvert=big-endian -O2 -fno-range-check -ffree-line-length-none
LFLAGS =

F90_VERSION = $(shell $(F90) --version | head -1)

# option to treat default real as real*8
R8 = -fdefault-real-8 -fdefault-double-8
EXTENDED_SOURCE = -ffixed-line-length-132

#
# Set the following to ensure that the beginning/end of records
# in sequential-access unformatted files have 4-byte markers
# (this does impose a 2 GB limit on record sizes, however)
#
#FFLAGS += -frecord-marker=4
#F90FLAGS += -frecord-marker=4

# check if ABI was specified explicitly
ifneq ($(ABI),)
FFLAGS += -m$(ABI)
F90FLAGS += -m$(ABI)
LFLAGS += -m$(ABI)
endif


# machine-specific options
ifeq ($(MACHINE),IRIX64)
FFLAGS += -mabi=64 
F90FLAGS += -mabi=64
LFLAGS += -mabi=64
endif

# uncomment next two lines for extensive debugging
# the following switch adds extra debugging
ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += -fbounds-check -fcheck-array-temporaries -ffpe-trap=invalid,zero,overflow -fbacktrace
F90FLAGS += -fbounds-check -fcheck-array-temporaries -ffpe-trap=invalid,zero,overflow -fbacktrace
FFLAGS += -finit-real=snan
F90FLAGS += -finit-real=snan
#LFLAGS += -lefence
endif
