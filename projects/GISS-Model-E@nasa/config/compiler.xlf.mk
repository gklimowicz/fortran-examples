
#F90 = xlf90_r
F90 = f95
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
# ibm compiler doesn't understand "-D" . Have to use "-WF,-D..."
CPP_FLAG_PREFIX = '-WF,'
CPPFLAGS = -DMACHINE_IBM -DCOMPILER_XLF
FFLAGS = -O2 -qfixed -qsuffix=cpp=f -qmaxmem=16384
F90FLAGS = -O2 -qfree -qsuffix=cpp=F90 -qmaxmem=16384
# one may need to add -bmaxstack:0x1000000 if rusns out of stack
LFLAGS = -O2 # -bmaxdata:0x10000000
# no guarantee that the following line gives correct info
#F90_VERSION = $(shell $(F90) -qversion 2>&1)

# redefine values for IBM workstation
ifeq ($(MACHINE),AIX)
F90 = xlf90_r
CMP_MOD = $(SCRIPTS_DIR)/compare_module_file.pl -compiler IBM-xlf90-on-AIX
F90_VERSION = $(shell what /usr/lpp/xlf/bin/xlfentry | tail -1)
endif
