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

## The files in this directory implement objects and methods
## as they appear in the Sindarin language

## We create a library which is still to be combined with auxiliary libs.
noinst_LTLIBRARIES = libtypes.la
check_LTLIBRARIES = libtypes_ut.la

libtypes_la_SOURCES = \
  $(TYPES_MODULES) \
  $(TYPES_SUBMODULES)

TYPES_MODULES = \
  particle_specifiers.f90 \
  analysis.f90 \
  pdg_arrays.f90 \
  jets.f90 \
  subevents.f90

TYPES_SUBMODULES = \
  particle_specifiers_sub.f90 \
  analysis_sub.f90 \
  pdg_arrays_sub.f90 \
  subevents_sub.f90

libtypes_ut_la_SOURCES = \
  particle_specifiers_uti.f90 particle_specifiers_ut.f90 \
  analysis_uti.f90 analysis_ut.f90 \
  pdg_arrays_uti.f90 pdg_arrays_ut.f90 \
  jets_uti.f90 jets_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = types.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ${TYPES_MODULES:.f90=.$(FCMOD)}

# Submodules must not be included here
libtypes_Modules = ${TYPES_MODULES:.f90=} ${libtypes_ut_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(libtypes_Modules); do \
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
  ../utilities/Modules \
  ../testing/Modules \
  ../system/Modules \
  ../combinatorics/Modules \
  ../parsing/Modules \
  ../physics/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(libtypes_la_SOURCES) $(libtypes_ut_la_SOURCES)
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
Makefile.depend: $(libtypes_la_SOURCES) $(libtypes_ut_la_SOURCES)
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

AM_FCFLAGS = -I../basics -I../utilities -I../testing -I../system -I../combinatorics -I../parsing -I../physics -I../qft -I../fastjet

########################################################################
# For the moment, the submodule dependencies will be hard-coded
particle_specifiers_sub.lo: particle_specifiers.lo
analysis_sub.lo: analysis.lo
pdg_arrays_sub.lo: pdg_arrays.lo
subevents_sub.lo: subevents.lo

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

## Dependencies across directories and packages, if not automatically generated
$(libtypes_la_OBJECTS): \
  ../fastjet/fastjet.$(FCMOD)

## (Re)create F90 sources from NOWEB source.
if NOWEB_AVAILABLE

PRELUDE = $(top_srcdir)/src/noweb-frame/whizard-prelude.nw
POSTLUDE = $(top_srcdir)/src/noweb-frame/whizard-postlude.nw

types.stamp: $(PRELUDE) $(srcdir)/types.nw $(POSTLUDE)
	@rm -f types.tmp
	@touch types.tmp
	for src in $(libtypes_la_SOURCES) $(libtypes_ut_la_SOURCES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	@mv -f types.tmp types.stamp

$(libtypes_la_SOURCES) $(libtypes_ut_la_SOURCES): types.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f types.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) types.stamp; \
	fi

endif


########################################################################
## Non-standard cleanup tasks
## Remove sources that can be recreated using NOWEB
if NOWEB_AVAILABLE
maintainer-clean-noweb:
	-rm -f *.f90 *.c
endif
.PHONY: maintainer-clean-noweb

## Remove those sources also if builddir and srcdir are different
if NOWEB_AVAILABLE
clean-noweb:
	test "$(srcdir)" != "." && rm -f *.f90 *.c || true
endif
.PHONY: clean-noweb

## Remove F90 module files
clean-local: clean-noweb
	-rm -f types.stamp types.tmp
	-rm -f *.$(FCMOD)
if FC_SUBMODULES
	-rm -f *.smod *.sub
endif

## Remove backup files
maintainer-clean-backup:
	-rm -f *~
.PHONY: maintainer-clean-backup

## Register additional clean targets
maintainer-clean-local: maintainer-clean-noweb maintainer-clean-backup
