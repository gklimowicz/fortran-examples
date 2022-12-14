
F90 = f90
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend -H
CMP_MOD = $(SCRIPTS_DIR)/compare_module_file.pl -compiler MIPSPRO-f90-on-IRIX64
#FFLAGS = -ftpp -macro_expand -O2 -64 -mips4 -OPT:reorg_comm=off -w2 -OPT:Olimit=6500 -ansi -woff124 -woff52 
FFLAGS = -cpp -Wp,-P -O2 -64 -mips4 -OPT:reorg_comm=off -w2 -OPT:Olimit=6500 -ansi -woff124 -woff52 -I.
F90FLAGS = -cpp -Wp,-P -O2 -64 -mips4 -OPT:reorg_comm=off -w2 -OPT:Olimit=6500 -ansi -woff124 -woff52 -freeform -I.
LFLAGS = -64 -O2 -mips4 -lfastm -OPT:reorg_common=OFF
ifeq ($(MP),YES)
FFLAGS += -mp
F90FLAGS += -mp
LFLAGS += -mp
endif
# suppress some linker warnings if no verbose output
ifeq ($(VERBOSE_OUTPUT),NO)
LFLAGS += -LD_MSG:OFF=84,85,15,134
endif
# the following switch adds extra debugging
ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += -DEBUG:div_check=3 -DEBUG:subscript_check=ON -DEBUG:trap_uninitialized=ON
LFLAGS += -DEBUG:conform_check=YES -DEBUG:div_check=3 -DEBUG:subscript_check=ON -DEBUG:trap_uninitialized=ON
F90FLAGS += -DEBUG:div_check=3 -DEBUG:subscript_check=ON -DEBUG:trap_uninitialized=ON
LFLAGSF += -DEBUG:conform_check=YES -DEBUG:div_check=3 -DEBUG:subscript_check=ON -DEBUG:trap_uninitialized=ON
endif
# not sure if the following will help the debugging ...
# FFLAGS += -DEBUG:verbose_runtime=ON
# LFLAGS += -DEBUG:verbose_runtime=ON
F90_VERSION = $(shell $(F90) -version 2>&1)
