
# if MPI*DIR are not defined try to get them from MPIDIR
ifneq ($(MPIDIR),)
  ifneq ($(wildcard $(MPIDIR)/include/mpi.h),)
    MPIINCLUDEDIR ?= $(MPIDIR)/include
  else
    ifneq ($(wildcard $(MPIDIR)/include/openmpi/mpi.h),)
      MPIINCLUDEDIR ?= $(MPIDIR)/include/openmpi
    endif
  endif
  ifeq ($(MPIINCLUDEDIR),)
    $(error MPI distribution not found. Check settings in ~/.modelErc)
  endif
  MPILIBDIR ?= $(MPIDIR)/lib
  MPIBINDIR  ?= $(MPIDIR)/bin
endif

# if MPI*DIR are not yet defined try to get them from a module
MPIINCLUDEDIR ?= $(MPI_INCLUDE)
MPILIBDIR ?= $(MPI_LIB)
MPIBINDIR  ?= $(MPI_BIN)

ifneq ($(MPIINCLUDEDIR),)
  FFLAGS += -I$(MPIINCLUDEDIR)
  F90FLAGS += -I$(MPIINCLUDEDIR)
  CPPFLAGS += -I$(MPIINCLUDEDIR)
endif
ifneq ($(MPILIBDIR),)
  LIBS += -L$(MPILIBDIR)
endif

# if MPIBINDIR not defined, assume that mpirun is in the PATH
ifneq ($(MPIBINDIR),)
  MPIRUN = $(MPIBINDIR)/mpirun
else
  MPIRUN = mpirun
endif


# try to work around memory leak
CPPFLAGS += -DMPITYPE_LOOKUP_HACK

VER := $(subst ., ,$(word 4,$(shell $(MPIRUN) --version 2>&1)))
VER_MAJOR := $(word 1,$(VER))
VER_MINOR := $(word 2,$(VER))
ifneq (,$(filter 7 8 9 10,$(VER_MINOR))$(filter 2 3 4,$(VER_MAJOR)))
LIBS += -lmpi_mpifh -lmpi
else
LIBS += -lmpi_f77 -lmpi
# -lmpi_cxx - this library may be needed for ESMF (?)
endif

ifneq ($(shell uname),Darwin)
  LIBS += -lrt
endif

