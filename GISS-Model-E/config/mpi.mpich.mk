

# if MPI*DIR are not defined try to get them from MPIDIR
ifneq ($(MPIDIR),)
  ifneq ($(wildcard $(MPIDIR)/include/mpi.h),)
    MPIINCLUDEDIR ?= $(MPIDIR)/include
  else
    ifneq ($(wildcard $(MPIDIR)/include/mpich2-x86_64/mpi.h),)
      MPIINCLUDEDIR ?= $(MPIDIR)/include/mpich2-x86_64
    endif
  endif
  ifeq ($(MPIINCLUDEDIR),)
    $(error MPI distribution not found. Check settings in ~/.modelErc)
  endif
  MPILIBDIR ?= $(MPIDIR)/lib
endif

# if MPI*DIR are not yet defined try to get them from a module
MPIINCLUDEDIR ?= $(MPI_INCLUDE)
MPILIBDIR ?= $(MPI_LIB)

ifneq ($(MPIINCLUDEDIR),)
  FFLAGS += -I$(MPIINCLUDEDIR)
  F90FLAGS += -I$(MPIINCLUDEDIR)
  CPPFLAGS += -I$(MPIINCLUDEDIR)
endif
ifneq ($(MPILIBDIR),)
  LIBS += -L$(MPILIBDIR)
endif

MPIRUN = mpirun

MPICH_VER := $(subst ., ,$(word 5,$(shell $(MPIRUN) -version 2>&1)))
MPICH_VER_MAJOR := $(word 1,$(MPICH_VER))
ifeq ($(MPICH_VER_MAJOR),3)
  LIBS += -lfmpich -lmpich -lopa -lmpl
else
  #LIBS += -ldl -lpmpich  -lmpich -lstdc++
  LIBS += -lpmpich -lmpich -lopa -lmpl -lpthread
endif

