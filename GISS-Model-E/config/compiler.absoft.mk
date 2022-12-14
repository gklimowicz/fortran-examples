
F90 = f90
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend -h
CMP_MOD = $(SCRIPTS_DIR)/compare_module_file.pl -compiler ABSOFT-f95-on-LINUX
CPPFLAGS += -DCONVERT_BIGENDIAN -DCOMPILER_Absoft
FFLAGS = -O2 -YEXT_NAMES=LCS -B108
F90FLAGS = -O2 -f free
LFLAGS = -lf90math -lV77 -lU77 -YEXT_NAMES=LCS -B108
# uncomment next two lines for extensive debugging
# the following switch adds extra debugging
ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += -g -trap=INVALID,DIVBYZERO,OVERFLOW -B111
LFLAGS += 
endif

# hack for module comparison hack
ifeq ($(MACHINE),Darwin)
CMP_MOD = $(SCRIPTS_DIR)/compare_module_file.pl -compiler ABSOFT-f95-on-DARWIN
endif

# use system cpp for preprocessing
EXTERNAL_CPP = YES

# absoft has non-standard flag for .mod files
I = p
