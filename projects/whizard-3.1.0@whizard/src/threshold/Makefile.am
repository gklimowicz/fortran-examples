## Makefile.am -- Makefile for model-parameter modules in WHIZARD
##
## Process this file with automake to produce Makefile.in
#
# Copyright (C) 1999-2022 by
#     Wolfgang Kilian <kilian@physik.uni-siegen.de>
#     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
#     Juergen Reuter <juergen.reuter@desy.de>
#     with contributions from
#     cf. main AUTHORS file
#
# WHIZARD is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# WHIZARD is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
########################################################################

## Build the objects required for the SM ttbar threshold model.

noinst_LTLIBRARIES = libthreshold.la
check_LTLIBRARIES = libthreshold_ut.la

## We need all in the distribution

EXTRA_DIST = interpolation.f90  nr_tools.f90  \
	toppik.f  toppik_axial.f  ttv_formfactors.f90

libthreshold_la_SOURCES = \
  $(THRESHOLD_MODULES) \
  $(THRESHOLD_SUBMODULES)

THRESHOLD_MODULES = \
  interpolation.f90 \
  nr_tools.f90 \
  toppik.f \
  toppik_axial.f \
  ttv_formfactors.f90

THRESHOLD_SUBMODULES = \
  interpolation_sub.f90 \
  ttv_formfactors_sub.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = threshold.nw

libthreshold_ut_la_SOURCES = \
  ttv_formfactors_uti.f90 ttv_formfactors_ut.f90

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  interpolation.$(FCMOD) \
  nr_tools.$(FCMOD) \
  ttv_formfactors.$(FCMOD) \
  nr.$(FCMOD) \
  nrtype.$(FCMOD) \
  nrutil.$(FCMOD) \
  hypgeo_info.$(FCMOD) \
  ode_path.$(FCMOD)

# Submodules must not be included here
# Dump module names into file Modules
libthreshold_Modules = \
	${THRESHOLD_MODULES:.f90=} \
	${libthreshold_ut_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(libthreshold_Modules); do \
          echo $$module >> $@.new; \
        done
	@if diff $@ $@.new -q >/dev/null; then \
          rm $@.new; \
	else \
          mv $@.new $@; echo "Modules updated"; \
        fi
BUILT_SOURCES = Modules

## Fortran module dependencies
# Get module lists from other directories
module_lists = \
  ../basics/Modules \
  ../testing/Modules \
  ../utilities/Modules \
  ../system/Modules \
  ../physics/Modules \
  ../combinatorics/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(libthreshold_la_SOURCES) \
  $(libthreshold_ut_la_SOURCES)
Module_dependencies.sed: $(module_lists)
	@rm -f $@
	echo 's/, *only:.*//' >> $@
	echo 's/, *&//' >> $@
	echo 's/, *.*=>.*//' >> $@
	echo 's/$$/.lo/' >> $@
	for list in $(module_lists); do \
		dir="`dirname $$list`"; \
		for mod in `cat $$list`; do \
			echo 's!: '$$mod'.lo$$!': $$dir/$$mod'.lo!' >> $@; \
		done \
	done

DISTCLEANFILES = Module_dependencies.sed

# The following line just says
#    include Makefile.depend
# but in a portable fashion (depending on automake's AM_MAKE_INCLUDE
@am__include@ @am__quote@Makefile.depend@am__quote@

Makefile.depend: Module_dependencies.sed
Makefile.depend: $(libthreshold_la_SOURCES) \
  $(libthreshold_ut_la_SOURCES)
	@rm -f $@
	for src in $^; do \
	  module="`basename $$src | sed 's/\.f[90][0358]//'`"; \
	  grep '^ *use ' $$src \
	    | grep -v '!NODEP!' \
	    | sed -e 's/^ *use */'$$module'.lo: /' \
	          -f Module_dependencies.sed; \
	done > $@

DISTCLEANFILES += Makefile.depend

nrtype.$(FCMOD) nr.$(FCMOD): nr_tools.$(FCMOD)

SUFFIXES = .lo .$(FCMOD)

# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
	@:

AM_FCFLAGS = -I../basics -I../testing -I../utilities -I../system -I../physics -I../combinatorics

########################################################################
# For the moment, the submodule dependencies will be hard-coded
interpolation_sub.lo: interpolation.lo
ttv_formfactors_sub.lo: ttv_formfactors.lo

########################################################################
## Default Fortran compiler options

## Profiling
if FC_USE_PROFILING
AM_FCFLAGS += $(FCFLAGS_PROFILING)
endif

## OpenMP
if FC_USE_OPENMP
AM_FCFLAGS += $(FCFLAGS_OPENMP)
endif

## MPI
if FC_USE_MPI
AM_FCFLAGS += $(FCFLAGS_MPI)
endif

########################################################################
## Non-standard targets and dependencies

## (Re)create F90 sources from NOWEB source.
if NOWEB_AVAILABLE

PRELUDE = $(top_srcdir)/src/noweb-frame/whizard-prelude.nw
POSTLUDE = $(top_srcdir)/src/noweb-frame/whizard-postlude.nw

threshold.stamp: $(PRELUDE) $(srcdir)/threshold.nw $(POSTLUDE)
	@rm -f threshold.tmp
	@touch threshold.tmp
	for src in \
          $(libthreshold_la_SOURCES) \
          $(libthreshold_ut_la_SOURCES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	@mv -f threshold.tmp threshold.stamp

$(libthreshold_la_SOURCES) $(libthreshold_ut_la_SOURCES): threshold.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
		rm -f threshold.stamp; \
		$(MAKE) $(AM_MAKEFLAGS) threshold.stamp; \
	fi

endif

########################################################################
## Explicit dependencies

interpolation.$(FCMOD): interpolation.lo
nr_tools.$(FCMOD): nr_tools.lo
ttv_formfactors.$(FCMOD): ttv_formfactors.lo
nr.$(FCMOD): nr_tools.lo
nrtype.$(FCMOD): nr_tools.lo
nrutil.$(FCMOD): nr_tools.lo
hypgeo_info.$(FCMOD): nr_tools.lo
ode_path.$(FCMOD): nr_tools.lo

########################################################################
## Non-standard cleanup tasks
## Remove sources that can be recreated using NOWEB
if NOWEB_AVAILABLE
maintainer-clean-noweb:
	-rm -f *.f90
endif
.PHONY: maintainer-clean-noweb

## Remove those sources also if builddir and srcdir are different
if NOWEB_AVAILABLE
clean-noweb:
	test "$(srcdir)" != "." && rm -f *.f *.f90 || true
endif
.PHONY: clean-noweb

## Remove F90 module files
clean-local: clean-noweb
	-rm -f *.stamp *.tmp
	-rm -f *.$(FCMOD)
if FC_SUBMODULES
	-rm -f *.smod *.sub
endif

## Remove backup files
maintainer-clean-local:
	-rm -f *~

## Register additional clean targets
maintainer-clean-local: maintainer-clean-noweb maintainer-clean-backup
