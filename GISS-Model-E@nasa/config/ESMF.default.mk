# these are the options needed to compile the code with ESMF library
# (MPI support should be included separately)


CPPFLAGS += -DUSE_ESMF

ESMF_COMM = $(MPIDISTR)
ifeq ($(ESMF_COMM),intel)
ESMF_COMM = intelmpi
endif

# hack to get compatible compiler name for ESMF directory
ESMF_COMPILER ?= $(COMPILER)
ifeq ($(COMPILER),IRIX64)
ESMF_COMPILER = default
endif

# the following variables specify location of ESMF library and includes
# they can be overwritten in ~/.modelErc file if necessary
ifeq ($(ESMF5_DIR),)
ESMFINCLUDEDIR ?= ${BASELIBDIR5}/include/esmf
ESMFLIBDIR ?= ${BASELIBDIR5}/lib
else
ESMFINCLUDEDIR ?= ${ESMF5_DIR}/mod/mod${ESMF_BOPT}/$(MACHINE).$(ESMF_COMPILER).64.$(ESMF_COMM).default
ESMFLIBDIR ?= ${ESMF5_DIR}/lib/lib${ESMF_BOPT}/$(MACHINE).$(ESMF_COMPILER).64.$(ESMF_COMM).default
endif

# the following tells make where to look for system mod files
VPATH += ${ESMFINCLUDEDIR}

FFLAGS += -I${ESMFINCLUDEDIR}
F90FLAGS += -I${ESMFINCLUDEDIR}
LIBS += ${ESMFLIBDIR}/libesmf.a -lstdc++

# if we don't have netcdf library add netcdf_stubs
ifndef NETCDFHOME
LIBS +=  -lnetcdf_stubs
endif
