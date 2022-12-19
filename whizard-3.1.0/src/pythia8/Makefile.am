## Makefile.am -- Makefile for WHIZARD
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

## The files in this directory implement the LHA user process interface for WHIZARD
## and the interface to Pythia8.

## The files in this directory end up in an auxiliary libtool library.
noinst_LTLIBRARIES = libwo_pythia8.la
check_LTLIBRARIES = libwo_pythia8_ut.la
if PYTHIA8_AVAILABLE
libwo_pythia8_la_SOURCES = \
	LHAWhizard_events.h \
	LHAWhizard.h \
	LHAWhizard.cpp \
	LHAWhizardWrap.cpp \
	Pythia8Wrap.h \
	Pythia8Wrap.cpp
libwo_pythia8_la_CPPFLAGS = $(PYTHIA8_CXXFLAGS)
else
libwo_pythia8_la_SOURCES = \
	LHAWhizard_dummy.f90 \
	Pythia8Wrap_dummy.f90
endif

libwo_pythia8_la_SOURCES += \
	whizard_lha.f90 \
	pythia8.f90

libwo_pythia8_ut_la_SOURCES = \
	whizard_lha_uti.f90 \
	whizard_lha_ut.f90 \
	pythia8_uti.f90 \
	pythia8_ut.f90

SRC_FROM_NOWEB = \
	whizard_lha.f90 \
	pythia8.f90 \
	whizard_lha_uti.f90 whizard_lha_ut.f90 \
	pythia8_uti.f90 pythia8_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = pythia8.nw


# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  pythia8.$(FCMOD) whizard_lha.$(FCMOD)

libwo_pythia8_Modules = ${libwo_pythia8_la_SOURCES:.f90=} ${libwo_pythia8_ut_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(libwo_pythia8_Modules); do \
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
	../events/Modules \
	../utilities/Modules \
	../particles/Modules \
	../physics/Modules \
	../qft/Modules \
  ../rng/Modules \
	../types/Modules \
	../testing/Modules \
	../system/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(libwo_pythia8_la_SOURCES) $(libwo_pythia8_ut_la_SOURCES)
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
		done ; \
	done

DISTCLEANFILES = Module_dependencies.sed

# The following line just says
#    include Makefile.depend
# but in a portable fashion (depending on automake's AM_MAKE_INCLUDE
@am__include@ @am__quote@Makefile.depend@am__quote@

Makefile.depend: Module_dependencies.sed
Makefile.depend: $(libwo_pythia8_la_SOURCES) ${libwo_pythia8_ut_la_SOURCES}
	@rm -f $@
	for src in $^; do \
		module="`basename $$src | sed 's/\.f[90][0358]//'`"; \
		grep '^ *use ' $$src \
		| grep -v '!NODEP!' \
		| sed -e 's/^ *use */'$$module'.lo: /' \
				-f Module_dependencies.sed; \
	done > $@

DISTCLEANFILES += Makefile.depend


# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
	@:
#	touch $@

AM_FCFLAGS = -I../basics -I../events -I../utilities -I../testing -I../system -I../physics -I../particles -I../qft -I../rng -I../types -I../combinatorics -I../fastjet -I../variables

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

FILTER =
if FC_USE_MPI
FILTER += -filter "sed 's/defn MPI:/defn/'"
endif

PRELUDE = $(top_srcdir)/src/noweb-frame/whizard-prelude.nw
POSTLUDE = $(top_srcdir)/src/noweb-frame/whizard-postlude.nw

Pythia8.stamp: $(PRELUDE) $(srcdir)/pythia8.nw $(POSTLUDE)
	@rm -f Pythia8.tmp
	@touch Pythia8.tmp
	for src in $(SRC_FROM_NOWEB); do \
		$(NOTANGLE) -R[[$$src]] $(FILTER) $^ | $(CPIF) $$src; \
		done
	@mv -f Pythia8.tmp Pythia8.stamp

$(SRC_FROM_NOWEB): Pythia8.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
		rm -f Pythia8.stamp; \
		$(MAKE) $(AM_MAKEFLAGS) Pythia8.stamp; \
	fi

endif


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
	test "$(srcdir)" != "." && rm -f *.f90 || true
endif
.PHONY: clean-noweb

## Remove F90 module files
clean-local: clean-noweb
	-rm -f Pythia8.stamp Pythia8.tmp
	-rm -f *.$(FCMOD)
if FC_SUBMODULES
	-rm -f *.smod
endif

## Remove backup files
maintainer-clean-backup:
	-rm -f *~
.PHONY: maintainer-clean-backup

## Register additional clean targets
maintainer-clean-local: maintainer-clean-noweb maintainer-clean-backup
