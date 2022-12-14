
F90 = lf95
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CMP_MOD = $(SCRIPTS_DIR)/compare_module_file.pl -compiler LAHEY-lf95-on-LINUX
CPPFLAGS += -DCONVERT_BIGENDIAN -DCOMPILER_Lahey
FFLAGS = -O -Cpp
LFLAGS = 
F90_VERSION = $(shell $(F90) --version | grep Release)
ifeq ($(MP),YES)
FFLAGS += --openmp
F90FLAGS += --openmp
LFLAGS += --openmp
endif
ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += --trap
LFLAGS += --trap
F90FLAGS += --trap
LFLAGSF += --trap
endif
