
F90 = f90
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend -H
FFLAGS = -cpp -O2 -mips4 -OPT:reorg_comm=off -w2 -OPT:Olimit=5745
F90FLAGS = -cpp -O2 -mips4 -OPT:reorg_comm=off -w2 -OPT:Olimit=5745 -freeform
LFLAGS = -O2 -mips4 -lfastm -mp -OPT:reorg_common=OFF -Wl,-woff,134 -Wl,-woff,15
F90_VERSION = $(shell $(F90) -version 2>&1)
