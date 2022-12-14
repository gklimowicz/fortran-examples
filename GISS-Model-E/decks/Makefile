.SUFFIXES:
.PHONY: force rundeck depend gcm setup clean vclean newstart exe cmp .R htmldoc setup_script setup_nocomp tests check_vars check_RUNSRC update_decks_repo

# main definitions
USER = $(shell whoami)
HOMEDIR = $(wildcard ~)
MODELERC ?= $(HOMEDIR)/.modelErc
#default number of processors for setup
NPROC = 1
NPES = 1
# configuration file with user-defined options:
ifeq ($(findstring $(MAKECMDGOALS),config configure),)
  ifeq ($(wildcard $(MODELERC)),)
    $(error $(MODELERC) not found. Was "make config" done for this account?)
  endif
  include $(MODELERC)
  ifeq ($(COMPILER),)
    $(error you need to specify COMPILER in your ~.modelErc)
  endif
endif

# modelE directory structure
## MODEL_E_ROOT ?= $(shell pwd)/..
MODEL_E_ROOT ?= $(abspath $(shell pwd)/..)
ifeq ($(BUILD_OUT_OF_SOURCE),YES)
  MODEL_E_BUILD_DIR = $(RUN)_obj
else
  MODEL_E_BUILD_DIR = $(MODEL_E_ROOT)
endif
MODEL_DIR = $(MODEL_E_BUILD_DIR)/model
AUX_DIR = $(MODEL_E_ROOT)/aux
SCRIPTS_DIR = $(MODEL_E_ROOT)/exec
DECKS_DIR = $(shell pwd)
BIN_DIR = $(DECKS_DIR)/$(RUN)_bin
CONFIG_DIR = $(MODEL_E_ROOT)/config
TESTS_DIR = $(MODEL_E_ROOT)/tests
TEMPLATES_DIR = $(MODEL_E_ROOT)/templates
HTMLDOC_DIR = $(DECKS_DIR)/$(RUN)_html
export HTMLDOC_DIR
RUN_COMMAND = $(SCRIPTS_DIR)/runE

RUNSRC  =
RUN     =
RUNDECKSRC = $(TEMPLATES_DIR)/$(RUNSRC).R
PRUNDECKSRC = $(TEMPLATES_DIR)/$(RUNSRC).py
RUNDECK = $(DECKS_DIR)/$(RUN).R

export COMPILER
export MODELERC
export MPI
export COMPONENTS
export PNETCDFHOME
export NETCDFHOME
# For unit tests
export F90FLAGS
export CPPFLAGS

# export for build process
export BUILD_OUT_OF_SOURCE
export MODEL_E_ROOT
export CONFIG_DIR
export SCRIPTS_DIR

ifeq ($(VERBOSE_OUTPUT),NO)
  MAKEFLAGS=-s
endif

help:
	$(MAKE) -C $(MODEL_DIR) help

sinclude $(CONFIG_DIR)/rules.mk
sinclude $(DECKS_DIR)/$(RUN).mk

ifort_version:
	@echo $(IFORT_RELEASE)

linkNonProduction:
	@if [ -e $(TEMPLATES_DIR)/nonProduction/$(RUNSRC).R ]; then \
	ln -sf $(TEMPLATES_DIR)/nonProduction/$(RUNSRC).R $(TEMPLATES_DIR); fi

rundeck: linkNonProduction check_RUNSRC
	@if [ -f $(DECKS_REPOSITORY)/$(RUN).R -a \
	-z "`find $(DECKS_REPOSITORY)/$(RUN).R -user $(USER)`" ]; then \
	echo 'name $(RUN) is already used by somebody else'; \
	echo 'please choose a unique one'; exit 1  ; fi
	@if [ -s $(DECKS_REPOSITORY)/$(RUN).R -a $(OVERWRITE) != YES ];\
	then echo '$(RUN).R already exists in repository'; \
	echo 'use OVERWRITE=YES to overwrite'; \
	exit 1; fi
	#$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) RUNSRC=$(RUNSRC)
	@if [ "$(RUN)" = "" ]; then \
	  echo 'You should specify run name on a command line (i.e. RUN=...)'; \
	  exit 1; fi
	@if [ -s $(PRUNDECKSRC) ]; then \
	  echo "Generating rundeck from $(PRUNDECKSRC)" ; \
	  python $(PRUNDECKSRC) > $(RUNDECK); \
	elif [ -s $(RUNDECKSRC) ]; then \
	  echo '---- Regenerating Rundeck $(RUNDECK) from default ($(RUNDECKSRC)) ----' ; \
	  sed 's/^'$(notdir $(RUNSRC))'\([:\. ]\)/'$(RUN)'\1/g' $(RUNDECKSRC) \
	  | $(SCRIPTS_DIR)/protect_ccp_options.pl -c | $(CPP) -I$(TEMPLATES_DIR) \
	  | $(SCRIPTS_DIR)/protect_ccp_options.pl -c -u \
	  | $(SCRIPTS_DIR)/adapt_to_local.pl \
	  > $(RUNDECK); \
	else \
	  echo "No source rundeck for $(RUNSRC)" ; exit 1; \
	fi
	@echo '---------             Rundeck  created               ---------' 
	@cp -f $(RUN).R $(DECKS_REPOSITORY)/$(RUN).R
	@chmod 644 $(DECKS_REPOSITORY)/$(RUN).R      #just in case ...
	@echo '### Information below describes your run. Do not delete! ###'\
	    >> $(DECKS_REPOSITORY)/$(RUN).R

depend: $(RUN).mk $(MODEL_DIR)/Makefile
	$(MAKE) check_sanity RUN=$(RUN)
	$(MAKE) -C $(MODEL_DIR) depend_all RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)

ifeq ($(BUILD_OUT_OF_SOURCE),YES)
$(MODEL_DIR)/Makefile: $(MODEL_E_ROOT)/model/Makefile
	if [ ! -d $(MODEL_E_BUILD_DIR) ] ; then  mkdir $(MODEL_E_BUILD_DIR) ; fi
	if [ ! -d $(MODEL_DIR) ] ; then  mkdir $(MODEL_DIR) ; fi
	cp $(MODEL_E_ROOT)/model/Makefile $(MODEL_DIR)/Makefile
endif


gcm $(BIN_DIR)/$(RUN).exe: $(RUN).mk $(BIN_DIR) depend check_vars
	$(MAKE) -C $(MODEL_DIR) gcm RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)
	-mv $(MODEL_DIR)/$(RUN).bin $(BIN_DIR)/$(RUN).exe
	@echo '----    Looks like compilation finished successfully     ----'

setup: gcm
	$(MAKE) setup_script RUN=$(RUN)
ifeq ($(MAKECMDGOALS),setup)
	@echo '--- The directory $(RUN) was created and prepared for the run ---'
	@echo '---                                                        ---'
	@echo '--- To start the first hour run from the initial conditions --'
	@echo '--- you may want to execute:                               ---'
	@echo '---                                                        ---'
	@echo '   ../exec/runE $(RUN) -cold-restart -np <NPROC>          '
	@echo '---                                                        ---'
	@echo '--- or use your favorite script which does the same        ---'
	@echo '--- You can use "make setup-run ... " to combine setup     ---'
	@echo '--- and first hour run.                                    ---'
	@echo
endif

tests:
ifeq ($(PFUNIT),)
	@echo '--- To run tests set PFUNIT environment variable.    ---'
else
	$(MAKE) -C $(TESTS_DIR) all RUN=$(RUN)
endif

setup_nocomp:
	@[ -d $(BIN_DIR) ] || mkdir $(BIN_DIR)
	@if [ ! -e $(BIN_DIR)/$(RUN).exe ]; then \
	  cp $(CMRUNDIR)/$(RUN)/$(RUN).exe $(BIN_DIR)/$(RUN).exe || \
	  { echo "No executable in run directory"; \
	  echo "You have to run 'gmake setup ...' first"; exit 1; } \
	fi
	$(MAKE) setup_script RUN=$(RUN)

setup_script:
	@echo '----- Saving Rundeck and other info to global repository -----'
	@if [ ! -d $(DECKS_REPOSITORY) ]; then \
	 echo '--- $(DECKS_REPOSITORY) does not exist. Rundeck not saved ---';\
	else \
	  if [ -f $(DECKS_REPOSITORY)/$(RUN).R -a \
	    -z "`find $(DECKS_REPOSITORY)/$(RUN).R -user $(USER)`" ]; then \
	    echo 'name $(RUN) is already used by somebody else'; \
	    echo 'please choose a unique one'; exit 1; \
	  else \
	    cp -f $(RUN).R $(DECKS_REPOSITORY)/$(RUN).R ;\
	    chmod 644 $(DECKS_REPOSITORY)/$(RUN).R  ;    \
	    echo '### Information below describes your run. Do not delete! ###'\
	    >> $(DECKS_REPOSITORY)/$(RUN).R ; \
	    date >> $(DECKS_REPOSITORY)/$(RUN).R ; \
	    echo "$(F90_VERSION)" >> $(DECKS_REPOSITORY)/$(RUN).R ;\
	    if [ -d $(MODEL_E_ROOT)/.git ]; then \
	      echo "Git reference: `git rev-parse HEAD 2>&1`" \
	      >> $(DECKS_REPOSITORY)/$(RUN).R; else \
	      echo 'Version: `cat ../model/version`, unknown revision' \
	      >> $(DECKS_REPOSITORY)/$(RUN).R ;\
	    fi ;\
	  fi ;\
	fi
	@echo '---------        Starting setup for $(RUN)        ----------'
	@echo '--------------------------------------------------------------'
ifeq ($(MPI),YES)
	@$(SETUP) $(SETUP_FLAGS) -mpi $(RUN)
else
	@$(SETUP) $(SETUP_FLAGS) $(RUN)
endif
	@echo '---------       Setup  finished successfully         ---------'


aux rvr csrvr CMPE002: depend $(BIN_DIR)
	$(MAKE) -C $(MODEL_DIR) gcmlib RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)
	$(MAKE) -C $(AUX_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)

auxinput: aux
	$(MAKE) input -C $(AUX_DIR) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)

auxqflux:
	$(MAKE) -C $(MODEL_DIR) gcmlib RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)
	$(MAKE) qflux -C $(AUX_DIR) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)

auxdeep: aux
	$(MAKE) odeep -C $(AUX_DIR) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)

ifeq ($(BUILD_OUT_OF_SOURCE),YES)
clean: check_RUN
endif

clean vclean clean_all:
	$(MAKE) -C $(MODEL_DIR) $@ $(OPTS_MAIN)
ifneq ($(PFUNIT),)
	$(MAKE) -C $(TESTS_DIR) clean
endif
	$(MAKE) -C $(AUX_DIR) $@ $(OPTS_MAIN)

exe: $(RUN)/$(RUN) gcm
	@if [ ! -s $(CMRUNDIR)/$(RUN)/I ] ; then \
	echo \
	  'directory $(CMRUNDIR)/$(RUN) or input files inside it are missing';\
	echo 'you should make "setup" first'; \
	exit 1; fi
	@mv -f $(BIN_DIR)/$(RUN).exe $(CMRUNDIR)/$(RUN)/ ; \
	if [ $$? -ne 0 ] ; then \
	  echo "Hmm. Couldn't write the executable to the run directory." ; \
	  echo "Maybe you don't have write permissions?" ; \
	  exit 1 ; fi
	@echo '--- Successfully updated the executable in the run directory ---'
	@if [ ! -d $(DECKS_REPOSITORY) ] ; then \
	  echo '--- $(DECKS_REPOSITORY) does not exist. Rundeck not saved ---'; \
	else \
	  $(MAKE) update_decks_repo RUN=$(RUN); fi

update_decks_repo:
	@if [ -n "`find $(DECKS_REPOSITORY)/$(RUN).R -user $(USER)`" ] ; then \
	echo 'Rundeck in repository is ok - will add info'; \
	else \
	if [ -f $(DECKS_REPOSITORY)/$(RUN).R ] ; then \
	  echo 'repository file for $(RUN) belongs to somebody else'; \
	  echo 'Info not written. Please check!'; exit 1 ; fi ; \
	echo 'repository file for $(RUN) is missing'; \
	echo 'Either you didn''t run setup yet or repository is corrupted'; \
	echo 'PLEASE CHECK! You may be loosing important data!'; echo ; \
	echo 'trying to create a new repository file'; \
	echo '--- no information on original rundeck ---' \
	  > $(DECKS_REPOSITORY)/$(RUN).R ; \
	cat $(RUN).R >> $(DECKS_REPOSITORY)/$(RUN).R ; \
	echo '### Information below describes your run.' \
	  >> $(DECKS_REPOSITORY)/$(RUN).R ; \
	if [ ! -f $(DECKS_REPOSITORY)/$(RUN).R ] ; then \
	  echo 'no luck...' ; exit 1 ; fi ; \
	echo 'new $(DECKS_REPOSITORY)/$(RUN).R created'; fi
	@echo 'adding corresponding info to rundeck repository'
	@perl -e "while(<>){ if(/^### Information/){exit;} print; }" \
	<$(DECKS_REPOSITORY)/$(RUN).R  >$(DECKS_REPOSITORY)/$(RUN).R.tmp
	@echo '### Executable was recompiled ' >> $(DECKS_REPOSITORY)/$(RUN).R
	@date >> $(DECKS_REPOSITORY)/$(RUN).R
	@echo "$(F90_VERSION)" >> $(DECKS_REPOSITORY)/$(RUN).R
	@if [ -d $(MODEL_E_ROOT)/.git ]; then \
	echo "Git reference: `git rev-parse HEAD 2>&1`" \
	>> $(DECKS_REPOSITORY)/$(RUN).R; else \
	echo 'Version: `cat ../model/version`, unknown revision' \
	>> $(DECKS_REPOSITORY)/$(RUN).R; fi
	@echo 'diff old.R new.R :' >> $(DECKS_REPOSITORY)/$(RUN).R
	@diff $(DECKS_REPOSITORY)/$(RUN).R.tmp $(RUN).R \
	>> $(DECKS_REPOSITORY)/$(RUN).R || echo 'rundeck changed'
	@rm -f $(DECKS_REPOSITORY)/$(RUN).R.tmp
	@echo done

setup-run: setup
	$(RUN_COMMAND) $(RUN) -cold-restart -np $(NPES) &

newstart:
	@if [ "$(RUN)" = "" ]; then \
	echo 'You should specify run name on a command line (i.e. RUN=...)'; \
	exit 1; fi
	rm -i $(CMRUNDIR)/$(RUN)/*

htmldoc: depend
	-rm -rf $(HTMLDOC_DIR)
	$(MAKE) -C $(MODEL_DIR) htmldoc_all RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)

config configure:
	@echo
	@if [ -z "$(ModelE_Support)" ]; then \
	echo "Please specify location for ModelE support directory with" ; \
	echo "    ModelE_Support=<path_to_support_directory>" ; exit 1 ; fi
	@echo "This will create a global configuration file for user `whoami`"
	@if [ -f $(MODELERC) -a "$(OVERWRITE)" != YES ]; then \
	echo '~/.modelErc already exists, use OVERWRITE=YES to overwrite'; \
	exit 1; fi
	@cp -f $(SCRIPTS_DIR)/modelErc $(MODELERC)
	@sed "s|<ModelE_Support>|$(ModelE_Support)|" < $(MODELERC) > \
	$(MODELERC).tmp; mv -f $(MODELERC).tmp $(MODELERC)
	@echo
	@echo '$(MODELERC) has been created.'
	@echo 'You may want to edit it manually to customize the options.'
	@echo
	@if [ -a $(ModelE_Support) ]; then \
	echo "$(ModelE_Support) is present - will not modify" ; \
	else \
	mkdir $(ModelE_Support) ; \
	mkdir $(ModelE_Support)/prod_decks ; \
	mkdir $(ModelE_Support)/prod_runs ; \
	mkdir $(ModelE_Support)/prod_input_files ; \
	mkdir $(ModelE_Support)/exec ; \
	mkdir $(ModelE_Support)/huge_space ; \
	echo '$(ModelE_Support) has been created.' ; fi
	@echo 'Make sure that all directories specified in your $(MODELERC)'
	@echo 'really do exist on your system.'
	@echo

src_cpp_dir:
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)

etags:
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) $(OPTS_MAIN)

$(RUN).mk: $(DECKS_DIR)/$(RUN).R $(SCRIPTS_DIR)/r_to_mk.pl
	@if [ -z "$(RUN)" ]; then echo "Please specify RUN name"; exit 1; fi
	$(SCRIPTS_DIR)/r_to_mk.pl $< > $@

$(BIN_DIR):
	mkdir $(BIN_DIR)

check_sanity:
	@if [ "$(COMPONENTS)"xx = xx ] ; then \
	  echo "Your rundeck contains no Components" ; \
	  echo "Old style rundeck ?  ... aborting ..." ; \
	  exit 1 ; fi

check_RUN:
	@if [ "$(RUN)"xx = xx ] ; then \
	  echo "Your have to specify rundeck name with RUN=..." ; \
	  echo "... aborting ..." ; \
	  exit 1 ; fi

check_RUNSRC:
	@if [ "$(RUNSRC)"xx = xx ] ; then \
	  echo "Your have to specify template name with RUNSRC=..." ; \
	  echo "For example: RUNSRC=E4F40" ; \
	  echo "... aborting ..." ; \
	  exit 1 ; fi

check_vars:
ifeq ($(ESMF),YES)
ifeq ($(FVCUBED),)
	@echo
	@echo '--- Rundeck $(RUN) does not work with ESMF=YES option. ---'
	@echo    
	@exit 1
endif
ifeq ($(FVCUBED),NO)
	@echo    
	@echo '--- Rundeck $(RUN) does not work with ESMF=YES option. ---'
	@echo
	@exit 1
endif
ifeq ($(ESMF5_DIR),)
ifeq ($(BASELIBDIR5),)
	@echo
	@echo "You have to set ESMF5_DIR or BASELIBDIR5 in your ~/.modelErc"
	@echo    	
	@exit 1
endif
endif
endif

