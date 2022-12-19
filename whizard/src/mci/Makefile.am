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

## The files in this directory implement quantum field theory concepts
## such as model representation and quantum numbers.

## We create a library which is still to be combined with auxiliary libs.
noinst_LTLIBRARIES = libmci.la
check_LTLIBRARIES = libmci_ut.la

COMMON_F90 = \
  mci_base.f90 \
  iterations.f90 \
  integration_results.f90 \
  mci_none.f90 \
  mci_midpoint.f90 \
  mci_vamp.f90 \
  mci_vamp2.f90 \
  dispatch_mci.f90
MPI_F90 = \
  mci_vamp2_sub.f90_mpi
SERIAL_F90 = \
  mci_vamp2_sub.f90_serial
COMMON_SUBMODULES = \
  mci_base_sub.f90 \
  iterations_sub.f90 \
  integration_results_sub.f90 \
  mci_none_sub.f90 \
  mci_midpoint_sub.f90 \
  mci_vamp_sub.f90 \
  dispatch_mci_sub.f90
MCI_MODULES = \
  $(COMMON_F90)
MCI_SUBMODULES = \
  $(COMMON_SUBMODULES) \
  mci_vamp2_sub.f90

EXTRA_DIST = \
  $(COMMON_F90) \
  $(COMMON_SUBMODULES) \
  $(SERIAL_F90) \
  $(MPI_F90)

nodist_libmci_la_SOURCES = \
  $(MCI_MODULES) \
  $(MCI_SUBMODULES)

DISTCLEANFILES = mci_vamp2_sub.f90

if FC_USE_MPI
mci_vamp2_sub.f90: mci_vamp2_sub.f90_mpi
	-cp -f $< $@
else
mci_vamp2_sub.f90: mci_vamp2_sub.f90_serial
	-cp -f $< $@
endif

libmci_ut_la_SOURCES = \
  mci_base_uti.f90 mci_base_ut.f90 \
  iterations_uti.f90 iterations_ut.f90 \
  mci_none_uti.f90 mci_none_ut.f90 \
  integration_results_uti.f90 integration_results_ut.f90 \
  mci_midpoint_uti.f90 mci_midpoint_ut.f90 \
  mci_vamp_uti.f90 mci_vamp_ut.f90 \
  mci_vamp2_uti.f90 mci_vamp2_ut.f90 \
  dispatch_mci_uti.f90 dispatch_mci_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = mci.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ${MCI_MODULES:.f90=.$(FCMOD)}

# Submodules must not be included here
libmci_Modules = $(MCI_MODULES:.f90=) $(libmci_ut_la_SOURCES:.f90=)
Modules: Makefile
	@for module in $(libmci_Modules); do \
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
  ../rng/Modules \
  ../phase_space/Modules \
  ../variables/Modules \
  ../vegas/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(nodist_libmci_la_SOURCES) $(libmci_ut_la_SOURCES)
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

DISTCLEANFILES += Module_dependencies.sed

# The following line just says
#    include Makefile.depend
# but in a portable fashion (depending on automake's AM_MAKE_INCLUDE
@am__include@ @am__quote@Makefile.depend@am__quote@

Makefile.depend: Module_dependencies.sed
Makefile.depend: $(nodist_libmci_la_SOURCES) $(libmci_ut_la_SOURCES)
	@rm -f $@
	for src in $^; do \
	  module="`basename $$src | sed 's/\.f[90][0358]//'`"; \
	  grep '^ *use ' $$src \
	    | grep -v '!NODEP!' \
	    | sed -e 's/^ *use */'$$module'.lo: /' \
	          -f Module_dependencies.sed; \
	done > $@

DISTCLEANFILES += Makefile.depend

## Dependencies across directories and packages, if not automatically generated
mci_vamp.lo: ../../vamp/src/vamp.$(FCMOD)

# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
	@:
#	touch $@

AM_FCFLAGS = -I../basics -I../utilities -I../testing -I../system -I../combinatorics -I../rng -I../physics -I../fastjet -I../qft -I../matrix_elements -I../types -I../particles -I../beams -I../phase_space -I../expr_base -I../variables -I../../vamp/src -I../vegas

########################################################################
mci_base_sub.lo: mci_base.lo
iterations_sub.lo: iterations.lo
integration_results_sub.lo: integration_results.lo
mci_none_sub.lo: mci_none.lo
mci_midpoint_sub.lo: mci_midpoint.lo
mci_vamp_sub.lo: mci_vamp.lo
mci_vamp2_sub.lo: mci_vamp2.lo
dispatch_mci_sub.lo: dispatch_mci.lo

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

# MPI
if FC_USE_MPI
AM_FCFLAGS += $(FCFLAGS_MPI)
endif

########################################################################
## Non-standard targets and dependencies

## (Re)create F90 sources from NOWEB source.
if NOWEB_AVAILABLE

FILTER = -filter "sed 's/defn MPI:/defn/'"

PRELUDE = $(top_srcdir)/src/noweb-frame/whizard-prelude.nw
POSTLUDE = $(top_srcdir)/src/noweb-frame/whizard-postlude.nw

mci.stamp: $(PRELUDE) $(srcdir)/mci.nw $(POSTLUDE)
	@rm -f mci.tmp
	@touch mci.tmp
	for src in $(COMMON_F90) $(libmci_ut_la_SOURCES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	for src in $(COMMON_SUBMODULES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	for src in $(SERIAL_F90:.f90_serial=.f90); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src'_serial'; \
        done
	for src in $(MPI_F90:.f90_mpi=.f90); do \
	  $(NOTANGLE) -R[[$$src]] $(FILTER) $^ | $(CPIF) $$src'_mpi'; \
        done
	@mv -f mci.tmp mci.stamp

$(COMMON_F90) $(COMMON_SUBMODULES) $(SERIAL_F90) $(MPI_F90) $(libmci_ut_la_SOURCES): mci.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f mci.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) mci.stamp; \
	fi

endif


########################################################################
## Non-standard cleanup tasks
## Remove sources that can be recreated using NOWEB
if NOWEB_AVAILABLE
maintainer-clean-noweb:
	-rm -f *.f90 *.f90_serial *.f90_mpi *.c
endif
.PHONY: maintainer-clean-noweb

## Remove those sources also if builddir and srcdir are different
if NOWEB_AVAILABLE
clean-noweb:
	test "$(srcdir)" != "." && rm -f *.f90 *.f90_serial *.f90_mpi *.c || true
endif
.PHONY: clean-noweb

## Remove F90 module files
clean-local: clean-noweb
	-rm -f mci.stamp mci.tmp
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
