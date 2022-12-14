
F90 = nagfor
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CPPFLAGS += -DCOMPILER_NAG
FFLAGS = -fpp -O0 -maxcontin=100 -kind=byte -dusty -f2003 -convert=BIG_ENDIAN -mismatch_all
F90FLAGS = -fpp -O2 -free -kind=byte  -f2003 -convert=BIG_ENDIAN -mismatch_all
LFLAGS =
EXTENDED_SOURCE = -132

# the following switch adds extra debugging
ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += -C=all -gline
F90FLAGS += -C=all -gline
#LFLAGS += -lefence
endif
