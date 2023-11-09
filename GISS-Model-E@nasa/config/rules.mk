#
# this file contains rules shared by Makefiles in "model" and "aux"
#

.PHONY:
.DELETE_ON_ERROR:

ifdef MOD_DIR
  VPATH += $(MOD_DIR)
endif

######  Some user customizable settings:   ########

# EXTRA_FFLAGS specifies some extra flags you want to pass
# to Fortarn compiler, like
# -g        - include debugging information
# -listing  - create listings (.L)
EXTRA_FFLAGS =

# EXTRA_LFLAGS specifies some extra flags you want to pass
# to linker. Currently needed as a hack to compile hybrid MPI/OpenMP
# code to pass "-openmp" to linker
EXTRA_LFLAGS =

# hack to force compilation errors to be written to ERR files rather than scren
ifeq ($(OUTPUT_TO_FILES),YES)
  COMP_OUTPUT = > $*.ERR 2>&1 || { r=$$? ; cat $*.ERR ; exit $$r ; }
  LINK_OUTPUT = > $(RUN).ERR 2>&1 || { r=$$? ; cat $(RUN).ERR ; exit $$r ; }
else
  COMP_OUTPUT =
  LINK_OUTPUT =
endif

# if -s specified enable some extra messages
ifeq ($(findstring s,$(MFLAGS)),s)
  MSG =
else
  MSG = > /dev/null
endif

#
# starting machine - specific options
#

NO_COMMAND = @echo "*****  This architecture is not supported "; \
             echo "*****  or compiler is not specified properly."; \
             echo "*****  You have COMPILER=$(COMPILER)" ; exit 1;
F90 = $(NO_COMMAND)
CC ?= cc
FMAKEDEP = $(NO_COMMAND)
CMP_MOD = cmp -s
SETUP = $(SCRIPTS_DIR)/setup_e.pl
CPP = $(NO_COMMAND)
LIBS =
INCS =
F90_VERSION = 'Unknown compiler version'
ECHO_FLAGS =
CPPFLAGS =
# the following is default fortran flag for path to include and .mod files
# it is redefined later for compilers with non-standard flags
I = I
# by default assume that fortran compiler can do cpp
EXTERNAL_CPP = NO
# assume that C compiler understands basic gcc flags
CFLAGS = -O2
# check if ABI was specified
ifneq ($(ABI),)
  CFLAGS += -m$(ABI)
endif
# define default name for m4
M4 = m4
# default runlib
RANLIB = ranlib
# default name for lib directory for current ABI
ifeq ($(ABI),64)
  LIBABI = lib64
else
  LIBABI = lib
endif

# RFLAGS returns rundeck options for the current object (i.e. for $*)
# it has effect only of OBJ_LIST_O is defined
RFLAGS = $(shell perl \
-e '$$_="$(OBJ_LIST_O)"; m/\b$* *\|([^|]*)\|/; print " $$1";' \
)


# check consistency of compilation flags and redefine some
# flags if necessary
ifeq ($(ESMF),YES)
  MPI = YES
endif

ifeq ($(FVCORE),YES)
  MPI = YES
endif

ifeq ($(FVCUBED),YES)
  ESMF = YES
endif

ifeq ($(COSP_SIM),YES)
   CPPFLAGS += -DCOSP_SIM
   FFLAGS += -$(I)COSP_SIM
   F90FLAGS += -$(I)COSP_SIM
endif

# hack to keep Intel8 name valid (only temporarily)
ifeq ($(COMPILER),Intel8)
  $(error please set "COMPILER=intel" in your ~/.modelErc)
endif

# include machine-specific options
MACHINE = $(shell uname)
include $(CONFIG_DIR)/machine.$(MACHINE).mk

#include compiler-specific options
ifneq ($(COMPILER),)
  include $(CONFIG_DIR)/compiler.$(COMPILER).mk
endif

### HACK !! - add source dir to CPPFLAGS
#ifneq ($(SRC_DIR),)
#  CPPFLAGS += -I$(SRC_DIR)
#endif

#check if m4 is present
ifeq ($(shell $(M4) --version),)
  $(error compatible m4 preprocessor was not found on your system)
endif

ifeq ($(MPI),YES)
  CPPFLAGS += -DUSE_MPI
endif

# include ESMF library if necessary (sets CPPFLAGS appropriately)
ifeq ($(ESMF),YES)
  include $(CONFIG_DIR)/ESMF.default.mk
endif

ifeq ($(FVCUBED),YES)
  ifndef FVCUBED_ROOT
     FVCUBED_ROOT = false
  endif

  # Cubed-sphere requires FVCORE and MPP enabled
  FVCORE=YES
  #MPP=YES but current FVcubed has its own MPP already

  # Testing options: ADIABATIC
  ADIABATIC=YES

  CUBED_SPHERE=YES

  FVINC = -I$(FVCUBED_ROOT)/$(MACHINE)/include
  FVINCS = $(FVINC) $(FVINC)/MAPL_Base $(FVINC)/MAPL_cfio $(FVINC)/FVdycoreCubed_GridComp
  INCS += $(FVINCS)
  FVINCx = $(FVCUBED_ROOT)/$(MACHINE)/include
  FVINCSx = $(FVINCx):$(FVINCx)/MAPL_Base:$(FVINCx)/FVdycoreCubed_GridComp
  ifdef SYSTEM_MOD_DIRS
    SYSTEM_MOD_DIRS = $(SYSTEM_MOD_DIRS):$(FVINCSx)
  else
    SYSTEM_MOD_DIRS = $(FVINCSx)
  endif
  LIBS += -L$(FVCUBED_ROOT)/$(MACHINE)/lib -lMAPL_cfio -lMAPL_Base -lFVdycoreCubed_GridComp -lfvdycore -lGMAO_mpeu
  # this extra -lesmf would not be needed if the ESMF stuff came after this section
  LIBS += $(ESMFLIBDIR)/libesmf.a
  ifdef NETCDFHOME
    NETCDFLIB ?= -L$(NETCDFHOME)/lib -lnetcdf
    LIBS += $(subst ",,$(NETCDFLIB))
    #"
    NETCDFINCLUDE ?= -I$(NETCDFHOME)/include
    FFLAGS += $(NETCDFINCLUDE)
    F90FLAGS += $(NETCDFINCLUDE)
    INCS += $(NETCDFINCLUDE)
  endif

endif

# If using Fortuna2-5 w/HDF5
#ifeq ($(FVCUBED),YES)
#  LIBS += -lhdf5_hl -lhdf5 -lz -lm -lmfhdf -ldf -lsz -ljpeg -lm  -lmfhdf -ldf  -lcurl -lrt -lm -lz -lm
#endif

ifeq ($(FVCORE),YES)
  ifndef FVCORE_ROOT
     FVCORE_ROOT = false
  endif
  CPPFLAGS += -DUSE_FVCORE
  ifneq ($(FVCUBED),YES)
    FVINC = -I$(FVCORE_ROOT)/$(MACHINE)/include
    CPPFLAGS += -DFVCUBED_SKIPPED_THIS -DCREATE_FV_RESTART
    INCS += $(FVINC) $(FVINC)/GEOS_Base $(FVINC)/GEOS_Shared $(FVINC)/GMAO_gfio_r8 $(FVINC)/GMAO_cfio_r8 $(FVINC)/GMAO_pilgrim $(FVINC)/FVdycore_GridComp  -I$(BASELIBDIR)/include
    LIBS += -L$(FVCORE_ROOT)/$(MACHINE)/lib  -lFVdycore_GridComp  -lGMAO_pilgrim -lGMAO_gfio_r8 -lGMAO_cfio_r8 -lGEOS_Shared -lGEOS_Base -L$(BASELIBDIR)/lib
    LIBS += -L${BASELIBDIR}/lib -lesmf
  endif
endif

ifeq ($(SKIP_FV),YES)
  CPPFLAGS+=-DSKIP_FV
endif

ifeq ($(CUBED_SPHERE),YES)
  CPPFLAGS += -DCUBED_SPHERE
  INCS += -I$(FFTW_ROOT)/include
  LIBS += -L$(FFTW_ROOT)/lib -lfftw3
endif

ifeq ($(MPP),YES)
  CPPFLAGS += -DUSE_MPP
  # if using MPP installation on /usr/local
  FFLAGS += -I$(MPPDIR)/include
  F90FLAGS += -I$(MPPDIR)/include
  LIBS += -L$(MPPDIR)/lib -lfms_mpp_shared
  # MPPDIR is the path of the GEOS5 installation
  # if using GFDL installation within GEOS5
  #FFLAGS += -I$(MPPDIR)/include/GFDL_fms
  #F90FLAGS += -I$(MPPDIR)/include/GFDL_fms
  #LIBS += -L$(MPPDIR)/lib -lGFDL_fms

  #LIBS += -lfmpi -lmpi
endif

ifdef PNETCDFHOME
  LIBS += -L$(PNETCDFHOME)/lib -lpnetcdf
  FFLAGS += -I$(PNETCDFHOME)/include
  INCS += -I$(PNETCDFHOME)/include
endif

ifneq ($(CUBED_SPHERE),YES)

ifdef NETCDFHOME
  ifneq ($(wildcard $(NETCDFHOME)/include/netcdf.inc),)
    NETCDFINCLUDEDIR ?= $(NETCDFHOME)/include
  else
    ifneq ($(wildcard $(NETCDFHOME)/include/netcdf-3/netcdf.inc),)
      NETCDFINCLUDEDIR ?= $(NETCDFHOME)/include/netcdf-3
    else
      $(error NetCDF include files not found)
    endif
  endif

  ifneq ($(wildcard $(NETCDFHOME)/$(LIBABI)/libnetcdf*),)
    NETCDFLIBDIR ?= $(NETCDFHOME)/$(LIBABI)
  else
    NETCDFLIBDIR ?= $(NETCDFHOME)/lib
  endif
endif

ifdef NETCDFINCLUDEDIR
  FFLAGS += -I$(NETCDFINCLUDEDIR)
  F90FLAGS += -I$(NETCDFINCLUDEDIR)
  INCS += -I$(NETCDFINCLUDEDIR)
endif

ifdef NETCDFLIBDIR
  LIBS += -L$(NETCDFLIBDIR) -lnetcdf
  ifeq ($(wildcard $(NETCDFLIBDIR)/libnetcdff.*),)
    LIBS += -L$(NETCDFLIBDIR) -lnetcdf
  else
    LIBS += -L$(NETCDFLIBDIR) -L/opt/local/lib -lnetcdff -lnetcdf
  endif
endif

endif


ifdef GLINT2INCLUDEDIR
  FFLAGS += -I$(GLINT2INCLUDEDIR)
  F90FLAGS += -I$(GLINT2INCLUDEDIR)
  INCS += -I$(GLINT2INCLUDEDIR)
endif

ifdef GLINT2LIBDIR
  LIBS += -L$(GLINT2LIBDIR) -lglint2 -L/opt/local/lib -lproj -lblitz -lCGAL -lmpfr -lgmp
endif



ifeq ($(ADIABATIC),YES)
  CPPFLAGS += -DADIABATIC
endif


ifeq ($(MPI),YES)
  ifeq ($(MPIDISTR),)
    # unknown distribution, just trying to add a mpi librarry ...
    LIBS += -lmpi
  else
    include $(CONFIG_DIR)/mpi.$(MPIDISTR).mk
  endif
endif

#!! hack to deal with Intel "source_include" bug
# basically has to assume that all include files are in model/shared directory
INCS += -I$(MODEL_E_ROOT)/model/shared


CPPFLAGS += $(INCS)

# path to the modules dir if present
ifdef MOD_DIR
  FFLAGS += -$(I)$(MOD_DIR)
  F90FLAGS += -$(I)$(MOD_DIR)
endif

ifdef INCLUDE_DIR
  CPPFLAGS += -I$(INCLUDE_DIR)
endif

ifeq ($(COMPARE_MODULES_HACK),NO)
CMP_MOD = cmp -s
endif

# add fortran flags into a single string
ifeq ($(EXTERNAL_CPP),YES)
FFLAGS_ALL =  $(FFLAGS) $(EXTRA_FFLAGS)
F90FLAGS_ALL = $(F90FLAGS) $(EXTRA_FFLAGS)
else
# hack to deal with some compilers (xlf) not understanding -D etc.
ifneq ($(CPP_FLAG_PREFIX),)
FFLAGS_ALL =  $(FFLAGS) $(EXTRA_FFLAGS) $(addprefix $(CPP_FLAG_PREFIX),$(CPPFLAGS))
F90FLAGS_ALL = $(F90FLAGS) $(EXTRA_FFLAGS) $(addprefix $(CPP_FLAG_PREFIX),$(CPPFLAGS))
else
FFLAGS_ALL =  $(FFLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS)
F90FLAGS_ALL = $(F90FLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS)
endif
endif

ifdef SYSTEM_MOD_DIRS
VPATH += $(subst :, ,$(SYSTEM_MOD_DIRS))
endif

#
# Pattern  rules
#

*.mod: .current_options

%.mod:
	@echo checking $@
	@if [ "$<empty" = "empty" ]; then \
	echo "No dependency for $@ : assuming it is a system module";\
	else \
	cmp $< $@ ; \
	if ! cmp -s $< $@ ; then \
	  echo "will copy $< to $@" ; \
	  cp $< $@ ; \
	fi ; \
	fi

%.smod:
	@if [ ! -f $@ ] ; then rm -f $< ; $(MAKE) $< RUN=$(RUN); fi


# Standard fortran
ifeq ($(EXTERNAL_CPP),YES)
%.o: %.f.cpp.f
else
%.o: %.f
endif
	@echo $(ECHO_FLAGS)  compiling `basename $<` ... $(MSG) \\c
	$(F90) -c -o $@ $(FFLAGS_ALL) $(RFLAGS) $< $(COMP_OUTPUT)
	@if [ -s $(DEPENDFILE) ] ; then \
	for i in \
	`perl -e 'while(<>){ if(/(\S+)\.mod: *(\w+\@$*\.smod)/){print " $$1";} }' $(DEPENDFILE)` ; \
	do \
	   cmp $$i.mod $$i\@$*\.smod > /dev/null 2>&1 ; if [ $$? -ne 0 ] ; then \
	    echo "$$i.mod updated" ; \
	    cp -f $$i.mod $$i\@$*\.smod ; cp $$i\@$*\.smod $$i.mod ; \
	   else \
	    echo "$$i.mod not changed - will skip recompilations" ; \
	    touch -r $$i\@$*\.smod $$i.mod ; \
	   fi ; \
	done ; \
	fi
	@if [ -s $*.ERR ] ; then echo $(MSG); else echo Done $(MSG); fi
ifdef COMP_OUTPUT
	@if [ -s $*.ERR ] ; then cat $*.ERR; else  rm -f $*.ERR; fi
endif

ifeq ($(EXTERNAL_CPP),YES)
%.o: %.F90.cpp.F90
else
%.o: %.F90
endif
	@echo $(ECHO_FLAGS)  compiling `basename $<` ... $(MSG) \\c
	$(F90) -c -o $@ $(F90FLAGS_ALL) $(RFLAGS) $< $(COMP_OUTPUT)
	@if [ -s $(DEPENDFILE) ] ; then \
	for i in \
	`perl -e 'while(<>){ if(/(\S+)\.mod: *(\w+\@$*\.smod)/){print " $$1";} }' $(DEPENDFILE)` ; \
	do \
	   cmp $$i.mod $$i\@$*\.smod > /dev/null 2>&1 ; if [ $$? -ne 0 ] ; then \
	    echo "$$i.mod updated" ; \
	    cp -f $$i.mod $$i\@$*\.smod ; touch $$i.mod ; \
	   else \
	    echo "$$i.mod not changed - will skip recompilations" ; \
	    touch -r $$i\@$*\.smod $$i.mod ; \
	   fi ; \
	done ; \
	fi
	@if [ -s $*.ERR ] ; then echo $(MSG); else echo Done $(MSG); fi
ifdef COMP_OUTPUT
	@if [ -s $*.ERR ] ; then cat $*.ERR; else  rm -f $*.ERR; fi
endif

# cpp preprocessing

%.f.cpp.f: %.f
	@echo $(ECHO_FLAGS)  preprocessing $< ... $(MSG) \\c
	$(CPP) $(CPPFLAGS) $*.f | sed -n '/^#pragma/!p' > $@

%.F90.cpp.F90: %.F90
	@echo $(ECHO_FLAGS)  preprocessing $< ... $(MSG) \\c
	$(CPP) $(CPPFLAGS) $*.F90 | sed -n '/^#pragma/!p' > $@

%.f.cpp: %.f
	@#echo preprocessing $<  $(MSG)
	$(CPP) $(CPPFLAGS) $< > $@

%.F90.cpp: %.F90
	 @#echo preprocessing $<  $(MSG)
	 $(CPP) $(CPPFLAGS) $< > $@

%.o: %.c
	$(CC) -c -O2 -m64 $<

%.f: %.m4f
	-rm -f $@
	$(M4) -I`dirname $<` $< > $@
	chmod -w $@

%.F90: %.m4F90
	-rm -f $@
	$(M4) -I`dirname $<` $< > $@
	chmod -w $@



# end of Pattern  rules

