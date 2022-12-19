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

## The files in this directory make up the shower routines

noinst_LTLIBRARIES = libshower.la
check_LTLIBRARIES = libshower_ut.la

libshower_la_SOURCES = \
  $(SHOWER_MODULES) \
  $(SHOWER_SUBMODULES)

SHOWER_MODULES = \
  pythia6_up.f \
  ktclus.f90 \
  shower_base.f90 \
  shower_partons.f90 \
  shower_core.f90 \
  shower_pythia6.f90 \
  shower_pythia8.f90

SHOWER_SUBMODULES = \
  shower_base_sub.f90 \
  shower_partons_sub.f90 \
  shower_core_sub.f90 \
  shower_pythia6_sub.f90 \
  shower_pythia8_sub.f90

libshower_ut_la_SOURCES = \
  shower_base_uti.f90 shower_base_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = shower.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ktclus.$(FCMOD) \
  shower_base.$(FCMOD) \
  shower_core.$(FCMOD) \
  shower_partons.$(FCMOD) \
  shower_pythia6.$(FCMOD) \
  shower_pythia8.$(FCMOD)

libshower_Modules = ${SHOWER_MODULES:.f90=} ${libshower_ut_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(libshower_Modules); do \
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
  ../rng/Modules \
  ../physics/Modules \
  ../qft/Modules \
  ../expr_base/Modules \
  ../particles/Modules \
  ../types/Modules \
  ../variables/Modules \
  ../model_features/Modules \
  ../muli/Modules \
  ../events/Modules \
  ../beams/Modules \
  ../tauola/Modules \
  ../pythia8/Modules

include_modules_bare = ${module_lists:/Modules=}
include_modules = ${include_modules_bare:../%=-I../%}

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(libshower_la_SOURCES) $(libshower_ut_la_SOURCES)
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
Makefile.depend: $(libshower_la_SOURCES) $(libshower_ut_la_SOURCES)
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

AM_FCFLAGS = $(include_modules) -I../../vamp/src -I../fastjet -I../lhapdf -I ../pdf_builtin -I../tauola -I../../pythia6 -I../pythia8
AM_FFLAGS =

########################################################################
# For the moment, the submodule dependencies will be hard-coded
shower_base_sub.lo: shower_base.lo
shower_partons_sub.lo: shower_partons.lo
shower_core_sub.lo: shower_core.lo
shower_pythia6_sub.lo: shower_pythia6.lo
shower_pythia8_sub.lo: shower_pythia8.lo

########################################################################
## Default Fortran compiler options

## Profiling
if FC_USE_PROFILING
AM_FFLAGS += $(FCFLAGS_PROFILING)
AM_FCFLAGS += $(FCFLAGS_PROFILING)
endif

## OpenMP
if FC_USE_OPENMP
AM_FFLAGS += $(FCFLAGS_OPENMP)
AM_FCFLAGS += $(FCFLAGS_OPENMP)
endif

## MPI
if FC_USE_MPI
AM_FFLAGS += $(FCFLAGS_MPI)
AM_FCFLAGS += $(FCFLAGS_MPI)
endif

########################################################################
## Non-standard targets and dependencies

## (Re)create F90 sources from NOWEB source.

if NOWEB_AVAILABLE

PRELUDE = $(top_srcdir)/src/noweb-frame/whizard-prelude.nw
POSTLUDE = $(top_srcdir)/src/noweb-frame/whizard-postlude.nw

shower.stamp: $(PRELUDE) $(srcdir)/shower.nw $(POSTLUDE)
	@rm -f shower.tmp
	@touch shower.tmp
	for src in $(libshower_la_SOURCES) $(libshower_ut_la_SOURCES); do \
		$(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
	done
	@mv -f shower.tmp shower.stamp

$(libshower_la_SOURCES) $(libshower_ut_la_SOURCES): shower.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
		rm -f shower.stamp; \
		$(MAKE) $(AM_MAKEFLAGS) shower.stamp; \
	fi

endif

########################################################################
## Non-standard cleanup tasks
## Remove sources that can be recreated using NOWEB
if NOWEB_AVAILABLE
maintainer-clean-noweb:
	-rm -f *.f90 *.f
endif
.PHONY: maintainer-clean-noweb

## Remove those sources also if builddir and srcdir are different
if NOWEB_AVAILABLE
clean-noweb:
	test "$(srcdir)" != "." && rm -f *.f90 *.f || true
endif
.PHONY: clean-noweb

## Remove F90 module files
clean-local: clean-noweb
	-rm -f shower.stamp shower.tmp
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
