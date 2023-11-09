
F90 = ifort
IFORT_RELEASE := $(shell ifort --version | perl -e \
  'while(<>){ if(/ifort.* (\d+\.\d+)/) { print "$$1"; } }')
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CMP_MOD = $(SCRIPTS_DIR)/compare_module_file.pl -compiler INTEL-ifort-9-0-on-LINUX
FFLAGS = -fpp -O2 -ftz         -convert big_endian 
F90FLAGS = -fpp -O2 -ftz        -convert big_endian -free 
LFLAGS = -O2 -ftz
CPPFLAGS += -DCOMPILER_Intel8 -DCONVERT_BIGENDIAN
F90_VERSION = $(shell $(F90) -v 2>&1)
ifeq ($(MP),YES)
FFLAGS += -openmp
F90FLAGS += -openmp
LFLAGS += -openmp
endif
R8 = -r8
EXTENDED_SOURCE = -extend_source

SUPPORTED_RELEASES = 14.0 15.0 16.0 19.1

ifneq ($(OVERWRITE),YES)
ifeq ($(findstring $(IFORT_RELEASE),$(SUPPORTED_RELEASES)),)
  $(error ifort version $(IFORT_RELEASE) is not supported by this code. \
          Use one of: $(SUPPORTED_RELEASES) . \
          If you insist on using an unsupported version, you can do it at your \
          own risk by appending "OVERWRITE=YES" to the compilation command )
endif
endif


# flags needed for particular releases

# default flags for latest releases (work for 12.*, 13.*):
FFLAGS_RELEASE = -assume protect_parens -fp-model strict -warn nousage

# if some releases require different flags enter them here
ifeq ($(IFORT_RELEASE),11.1)
FFLAGS_RELEASE = -assume protect_parens -fp-model strict -fp-speculationoff
endif

ifeq ($(IFORT_RELEASE),11.0)
FFLAGS_RELEASE = -assume protect_parens -fp-model strict -fp-speculationoff
endif

ifeq ($(IFORT_RELEASE),10.1)
FFLAGS_RELEASE = -assume protect_parens -fp-model strict -diag-disable vec  -fp-speculationoff
endif

ifeq ($(IFORT_RELEASE),9.1)
FFLAGS_RELEASE =
endif

FFLAGS += $(FFLAGS_RELEASE)
F90FLAGS += $(FFLAGS_RELEASE)

ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += -CB -fpe0 -check uninit -ftrapuv -traceback
LFLAGS += -CB -fpe0 -check uninit -ftrapuv -traceback
F90FLAGS += -CB -fpe0 -check uninit -ftrapuv -traceback
LFLAGSF += -CB -fpe0 -check uninit -ftrapuv -traceback
endif
