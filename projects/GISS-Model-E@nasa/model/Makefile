################# NEW Makefile Format ###################################

.SUFFIXES:
.PHONY: FORCE

AVAILABLE_COMPONENTS = profiler MPI_Support Ent giss_LSM shared solvers dd2d tracers CS_Support

export MODEL_DIR := $(shell pwd)
export MOD_DIR := $(MODEL_DIR)/mod
export INCLUDE_DIR := $(MODEL_DIR)/include

export RUN_H = $(INCLUDE_DIR)/rundeck_opts.h

# export global options from OPTS_MAIN line here
export FVCUBED
export FVCORE

sinclude $(DECKS_DIR)/$(RUN).mk

all: $(COMPONENTS:=_dir)
	echo $^

define PROGRAM_template
$(1)_dir: $$($(1)_OBJS) $$($(1)_LIBS:%=-l%)
	@echo
	@echo "===> building component $(1)"
	$(MAKE) -C $(1) $(OPTS_$(notdir $(1)))
#	mv -f .liblist .liblist_tmp
	if [ -f $(1)/liblist ] ; then \
	  cat $(1)/liblist | sed 's,\blib,$(1)/lib,g' >> .liblist ; \
	else \
	  echo $(1)/lib$(notdir $(1)).a >> .liblist ; \
	fi
#	cat .liblist_tmp >> .liblist
#	rm -f .liblist_tmp
	@echo "===> component $(1) ok"
	@echo

$(1)_dep: $$($(1)_OBJS) $$($(1)_LIBS:%=-l%) $(RUN_H)
ifeq ($(BUILD_OUT_OF_SOURCE),YES)
	if [ ! -d $(1) ] ; then  mkdir $(1) ; fi
	cp $(MODEL_E_ROOT)/model/$(1)/Makefile $(1)
endif
	$(MAKE) -C $(1) depend $(OPTS_$(notdir $(1)))

$(1)_htmldoc:
	$(MAKE) -C $(1) htmldoc $(OPTS_$(notdir $(1)))

endef

$(foreach prog,$(COMPONENTS),$(eval $(call PROGRAM_template,$(prog))))


DEPENDFILE_SUBDIRS = .depend_subdirs
ifneq ($(MAKECMDGOALS),$(DEPENDFILE_SUBDIRS))
sinclude $(DEPENDFILE_SUBDIRS)
endif

depend_all: dir_structure $(COMPONENTS:=_dep)
	$(SCRIPTS_DIR)/comp_mkdep.pl $(COMPONENTS)
	+$(MAKE) depend

dir_structure:
	if [ ! -d $(INCLUDE_DIR) ] ; then  mkdir $(INCLUDE_DIR) ; fi
	if [ ! -d $(MOD_DIR) ] ; then  mkdir $(MOD_DIR) ; fi

TMPOBJSTRING = $(addsuffix .*,$(strip $(OBJ_LIST)))
ifeq ($(BUILD_OUT_OF_SOURCE),YES)
  TMPOBJSTRING :=  $(addprefix $(MODEL_E_ROOT)/model/,$(TMPOBJSTRING))
endif
FFFLIST = $(notdir $(wildcard $(TMPOBJSTRING)))
FSRCS = $(filter %.f, $(FFFLIST))
F90SRCS = $(filter %.F90, $(FFFLIST))

include $(CONFIG_DIR)/base.mk
#sinclude $(DEPENDFILE)
include $(CONFIG_DIR)/rules.mk

# sanity check - make sure that sources are present for all listed objects
MISSING_SRCS = $(filter-out $(FSRCS:.f=) $(F90SRCS:.F90=), $(strip $(OBJ_LIST)))
ifneq ($(MISSING_SRCS),)
  $(error no sources for requested objects: $(MISSING_SRCS))
endif

#COMPLIBS = $(patsubst %, %/lib.a, $(COMPONENTS))

#$(OBJS): $(COMPONENTS:=_dir)

COMPLIBS = $(shell perl -e 'print reverse <>;' < .liblist)

do_components: $(COMPONENTS:=_dir)

main.o: main.F90

do_main $(RUN).bin: main.o  $(OBJS) #  $(COMPONENTS:=_dir)
	-rm -f libmodel.a
	ar rcs libmodel.a $(OBJS)
	@echo "===> linking"
	$(F90) main.o $(LFLAGS) $(EXTRA_LFLAGS) $(OBJS) $(F90OBJS) $(ESMF_OBJS) \
	  $(COMPLIBS) $(LIBS) -o $(RUN).bin  $(LINK_OUTPUT)
	@echo "===> linking ok"
	@echo

main gcm: $(MOD_DIR)
	-rm -f .liblist
	touch .liblist
	+$(MAKE) do_components
	+$(MAKE) do_main

echo_vars:
	@echo CPP_OPTIONS = $(CPP_OPTIONS)
	@echo OBJ_LIST = $(OBJ_LIST)
	@echo COMPONENTS = $(COMPONENTS)
	@echo INPUT_FILES = $(INPUT_FILES)
	@echo RUN_PARAMETERS = $(RUN_PARAMETERS)
	@echo INPUTZ = $(INPUTZ)
	+$(MAKE) main

# the following line is for compatibility with old interface only
clean_all vclean: clean

# "clean" will clean components and then will do "clean" in main directory
clean: clean_components

clean_components:
	-rm -f $(RUN_H) $(MOD_DIR)/*.mod .liblist
	-for i in $(AVAILABLE_COMPONENTS) ; do \
	  if [ -d $$i ] ; then $(MAKE) -C $$i clean ; fi ; done

$(RUN_H): $(DECKS_DIR)/$(RUN).mk FORCE
	echo "#ifndef _RUNDECK_OPTS_H" > rundeck_opts.tmp ; \
	echo "#define _RUNDECK_OPTS_H" >> rundeck_opts.tmp ; \
	perl -e '$$_=qq{$(CPP_OPTIONS)}; s/ *\#/\n\#/g; print "$$_\n";' \
	 >> rundeck_opts.tmp ; \
	echo "#endif" >> rundeck_opts.tmp ; \
	if ! cmp -s rundeck_opts.tmp $(RUN_H) ; then \
	  mv rundeck_opts.tmp $(RUN_H) ; \
	else \
	  rm rundeck_opts.tmp ; \
	fi

$(DEPENDFILE): $(RUN_H)

$(LIB):

gcmlib: $(MOD_DIR)
	-rm -f .liblist
	touch .liblist
	+$(MAKE) do_components
	+$(MAKE) $(LIB)
#	mv -f .liblist .liblist_tmp
	echo $(LIB) >> .liblist
#	cat .liblist_tmp >> .liblist
#	rm -f .liblist_tmp

$(MOD_DIR):
	mkdir $(MOD_DIR)

# html documentation
htmldoc_all: $(COMPONENTS:=_htmldoc)
	+$(MAKE) htmldoc


