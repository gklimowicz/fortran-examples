# 
# this Makefile compiles programs in "aux" directory
#

.SUFFIXES:

.PHONY: clean vclean FORCE

MODEL_E_ROOT = ..
MODEL_DIR = $(MODEL_E_ROOT)/model
DECKS_DIR = $(MODEL_E_ROOT)/decks
SCRIPTS_DIR = $(MODEL_E_ROOT)/exec
ESMF_Interface=$(MODEL_DIR)/ESMF_Interface
CONFIG_DIR = $(MODEL_E_ROOT)/config
MOD_DIR = $(MODEL_DIR)/mod
INCLUDE_DIR = $(MODEL_DIR)/include

RUN     =
RUNDECK = $(DECKS_DIR)/$(RUN).R

OBJ_DIR = .
BIN_DIR = $(DECKS_DIR)/$(RUN)_bin

GCMLIB = $(addprefix $(MODEL_DIR)/,$(shell perl -e 'print reverse <>;' < $(MODEL_DIR)/.liblist))

#DEPENDFILE = $(OBJ_DIR)/.depend.$(RUN)

DEPENDFILE = NO_DEPENDFILE

# getting options from configuration files
# default modelE configuration file:
include $(SCRIPTS_DIR)/modelErc
# configuration file with user-defined options:
MODELERC ?= $(HOMEDIR)/.modelErc
sinclude $(MODELERC)


# list of targets to be built by default
all: $(BIN_DIR)/CMPE002 $(BIN_DIR)/CMPE002P $(BIN_DIR)/qc \
     $(BIN_DIR)/rsf2gic $(BIN_DIR)/ext_gic $(BIN_DIR)/RMS status
CMPE001:  $(BIN_DIR)/CMPE001
CMPE002:  $(BIN_DIR)/CMPE002
qflux: $(BIN_DIR)/vertflux.exe $(BIN_DIR)/ocnIC.exe $(BIN_DIR)/otspec.E001.exe
# $(BIN_DIR)/rsfE001.addO.exe $(BIN_DIR)/conv_rsf $(BIN_DIR)/conv_strat_rsf $(BIN_DIR)/conv_gic !! all obsolete
odeep: $(BIN_DIR)/mkdeep.exe
rvr: $(BIN_DIR)/RDconvert 
csrvr: RDijk2llEM
input: inputll2cs
regrid: ll2cs ncll2cs ncGinoux

# include file with implicit rules
include $(CONFIG_DIR)/rules.mk

# look for .mod files in MODEL_DIR
ifeq ($(COMPILER),Absoft)
FFLAGS += -p $(MOD_DIR) -I$(MODEL_DIR)
F90FLAGS += -p $(MOD_DIR) -I$(MODEL_DIR)
endif

ifeq ($(COMPILER),Intel)
FFLAGS += -convert big_endian
F90FLAGS += -convert big_endian
EXTRA_LFLAGS += -convert big_endian
else
  ifeq ($(COMPILER),gfortran)
  FFLAGS += -fconvert=big-endian
  F90FLAGS += -fconvert=big-endian
#  EXTRA_LFLAGS += -fconvert=big_endian
  endif
endif


#
#  targets
#

aux: all 

FORCE:

#$(BIN_DIR)/CMPE002: CMPE002.F90 $(GCMLIB)
#	$(F90) $(F90FLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

CMPE002.o: $(GCMLIB)

$(BIN_DIR)/CMPE002: CMPE002.o $(GCMLIB)
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

CMPE002P.o: $(GCMLIB)

$(BIN_DIR)/CMPE002P: CMPE002P.o $(GCMLIB)
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)
    
CMPE001.o: $(GCMLIB)

$(BIN_DIR)/CMPE001: CMPE001.o $(GCMLIB)
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

RMS.o: $(GCMLIB)

$(BIN_DIR)/RMS: RMS.o $(GCMLIB)
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

conv_rsf.o: $(GCMLIB)

$(BIN_DIR)/conv_rsf: conv_rsf.o $(GCMLIB)
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

conv_strat_rsf.o: $(GCMLIB)

$(BIN_DIR)/conv_strat_rsf: conv_strat_rsf.o $(GCMLIB)
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

conv_gic.o: $(GCMLIB)

$(BIN_DIR)/conv_gic: conv_gic.o $(GCMLIB)
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

rsf2gic.o: $(GCMLIB)

$(BIN_DIR)/rsf2gic: rsf2gic.o $(GCMLIB)
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

ext_gic.o: $(GCMLIB)

$(BIN_DIR)/ext_gic: ext_gic.o $(GCMLIB)
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

qc.o: $(GCMLIB)

$(BIN_DIR)/qc: qc.o $(GCMLIB)
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

vertflux.E001.o: $(GCMLIB) MAPS.o

$(BIN_DIR)/vertflux.exe: vertflux.E001.o MAPS.o $(GCMLIB)
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

ocnIC.E001.o: $(GCMLIB) MAPS.o

$(BIN_DIR)/ocnIC.exe: ocnIC.E001.o MAPS.o $(GCMLIB)
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

otspec.E001.o: $(GCMLIB) MAPS.o

$(BIN_DIR)/otspec.E001.exe: otspec.E001.o MAPS.o $(GCMLIB)
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

mkdeep.o: $(GCMLIB)

$(BIN_DIR)/mkdeep.exe: mkdeep.o $(GCMLIB)
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

RDconvert.o: $(GCMLIB)

$(BIN_DIR)/RDconvert: RDconvert.o $(GCMLIB)
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

RDijk2ll: RDijk2ll.o
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

RDijk2llEM: RDijk2llEM.o
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

inputll2cs: regrid.o ncio.o inputLL2CS.o 
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

ll2cs: regrid.o ll2cs.o  
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS) 

ncll2cs: regrid.o ncio.o ncll2cs.o
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS) 

ncGinoux: ncio.o ncGinoux.o
	$(F90) $(LFLAGS) $(EXTRA_LFLAGS) -o $@ $^ $(LIBS)

status:
	echo "AUX build successful"

clean: clean_all
	rm -f *.LST *~ *.ERR *.CHK *.L rundeck_opts.h

clean_all vclean:
	rm -f *.o *.mod *.smod .depend* .current_options

