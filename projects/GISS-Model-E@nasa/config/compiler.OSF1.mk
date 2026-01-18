
F90 = f90
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CMP_MOD = $(SCRIPTS_DIR)/compare_module_file.pl -compiler DEC-f90-on-OSF1
# CPPFLAGS = -DCONVERT_BIGENDIAN -DMACHINE_DEC
# FFLAGS = -O2 -cpp
# LFLAGS = -O2
FFLAGS = -O2 -cpp -convert big_endian
F90FLAGS = -O2 -cpp -convert big_endian -free
LFLAGS = -O2 -convert big_endian
ifeq ($(MP),YES)
FFLAGS += -omp
F90FLAGS += -omp
LFLAGS += -omp
endif
# the following switch adds extra debugging
ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += -check bounds -check overflow -fpe0
LFLAGS += -check bounds -check overflow -fpe0
F90FLAGS += -check bounds -check overflow -fpe0
LFLAGSF += -check bounds -check overflow -fpe0
endif
F90_VERSION = $(shell $(F90) -version 2>&1)
