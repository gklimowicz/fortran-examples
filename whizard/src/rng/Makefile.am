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

## The files in this directory implement standard algorithms for WHIZARD

## We create a library which is still to be combined with auxiliary libs.
noinst_LTLIBRARIES = librng.la
check_LTLIBRARIES = librng_ut.la

librng_la_SOURCES = \
  $(RNG_MODULES) \
  $(RNG_SUBMODULES)

RNG_MODULES = \
  rng_base.f90 \
  selectors.f90 \
  rng_tao.f90 \
  rng_stream.f90 \
  dispatch_rng.f90

RNG_SUBMODULES = \
  rng_base_sub.f90 \
  selectors_sub.f90 \
  rng_tao_sub.f90 \
  rng_stream_sub.f90 \
  dispatch_rng_sub.f90

librng_ut_la_SOURCES = \
  rng_base_uti.f90 rng_base_ut.f90 \
  selectors_uti.f90 selectors_ut.f90 \
  rng_tao_uti.f90 rng_tao_ut.f90 \
  rng_stream_uti.f90 rng_stream_ut.f90 \
  dispatch_rng_uti.f90 dispatch_rng_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = rng.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ${RNG_MODULES:.f90=.$(FCMOD)}

librng_Modules = ${RNG_MODULES:.f90=} ${librng_ut_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(librng_Modules); do \
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
  ../variables/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(librng_la_SOURCES) $(librng_ut_la_SOURCES)
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
Makefile.depend: $(librng_la_SOURCES) $(librng_ut_la_SOURCES)
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
rng_tao.lo: ../../vamp/src/tao_random_numbers.$(FCMOD)

# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
	@:
#	touch $@

AM_FCFLAGS = -I../basics -I../utilities -I../testing -I../system -I../../vamp/src -I../variables -I../expr_base -I../parsing -I../types -I../fastjet -I../physics -I../combinatorics

########################################################################
# For the moment, the submodule dependencies will be hard-coded
rng_base_sub.lo: rng_base.lo
selectors_sub.lo: selectors.lo
rng_tao_sub.lo: rng_tao.lo
rng_stream_sub.lo: rng_stream.lo
dispatch_rng_sub.lo: dispatch_rng.lo

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

rng.stamp: $(PRELUDE) $(srcdir)/rng.nw $(POSTLUDE)
	@rm -f rng.tmp
	@touch rng.tmp
	for src in $(librng_la_SOURCES) $(librng_ut_la_SOURCES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	@mv -f rng.tmp rng.stamp

$(librng_la_SOURCES) $(librng_ut_la_SOURCES): rng.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f rng.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) rng.stamp; \
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
	-rm -f rng.stamp rng.tmp
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
