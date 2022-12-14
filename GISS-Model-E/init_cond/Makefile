# 
# this Makefile compiles programs in "init_cond" directory
#

.SUFFIXES:

.PHONY: clean vclean


MODEL_E_ROOT = ..
MODEL_DIR = $(MODEL_E_ROOT)/model
DECKS_DIR = $(MODEL_E_ROOT)/decks
SCRIPTS_DIR = $(MODEL_E_ROOT)/exec
CONFIG_DIR = $(MODEL_E_ROOT)/config

# hack to prevent make from hanging
DEPENDFILE=FILE_DOES_NOT_EXIST

# getting options from configuration files
# default modelE configuration file:
include $(SCRIPTS_DIR)/modelErc
# configuration file with user-defined options:
sinclude ~/.modelErc         

# include file with implicit rules
sinclude $(CONFIG_DIR)/rules.mk


RES_FILE = UNDEFINED_RES_FILE
#
#  targets
#

# make sure resulution file is compiled first
AIC.D771201.o: $(RES_FILE).o

mkAIC.bin: $(RES_FILE).o AIC.D771201.o HNTRPS.o
	$(F90) $(FFLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

OIC.WOA98.exe: OIC.WOA98.o
	$(F90) $(FFLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

OPF.exe: OPF.o
	$(F90) $(FFLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

clean vclean:
	rm -f *.o *.mod *.vo *.inc .depend* *.sig *.sig1

UNDEFINED_RES_FILE.o:
	@echo you have to specify resolution file with RES_FILE=RES_XXX
	@echo use any of:
	@for i in `ls $(MODEL_DIR)/RES_*.f` ; do basename $${i%.*} ; done
	exit 1

# if file not present in current dir look at $(MODEL_DIR)
%.o: $(MODEL_DIR)/%.f
	$(F90) -c $(FFLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS) $(RFLAGS) $^ -o $@
