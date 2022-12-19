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
noinst_LTLIBRARIES = libme_methods.la
check_LTLIBRARIES = libme_methods_ut.la

libme_methods_la_SOURCES = \
  $(ME_METHODS_MODULES) \
  $(ME_METHODS_SUBMODULES)

ME_METHODS_MODULES = \
  prc_core.f90 \
  prc_test_core.f90 \
  prc_template_me.f90 \
  prc_omega.f90 \
  prc_external.f90 \
  prc_threshold.f90

ME_METHODS_SUBMODULES = \
  prc_core_sub.f90 \
  prc_test_core_sub.f90 \
  prc_template_me_sub.f90 \
  prc_omega_sub.f90 \
  prc_external_sub.f90 \
  prc_threshold_sub.f90

libme_methods_ut_la_SOURCES = \
  prc_template_me_uti.f90 prc_template_me_ut.f90 \
  prc_omega_uti.f90 prc_omega_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = me_methods.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ${ME_METHODS_MODULES:.f90=.$(FCMOD)}

libme_methods_Modules = \
  ${ME_METHODS_MODULES:.f90=}\
  ${libme_methods_ut_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(libme_methods_Modules); do \
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
  ../expr_base/Modules \
  ../physics/Modules \
  ../qft/Modules \
  ../types/Modules \
  ../matrix_elements/Modules \
  ../beams/Modules \
  ../variables/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(libme_methods_la_SOURCES) $(libme_methods_ut_la_SOURCES)
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
Makefile.depend: $(libme_methods_la_SOURCES) $(libme_methods_ut_la_SOURCES)
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

AM_FCFLAGS = -I../basics -I../utilities -I../testing -I../system -I../combinatorics -I../parsing -I../physics -I../qft -I../expr_base -I../types -I../matrix_elements -I../particles -I../beams -I ../variables -I../fastjet -I../pdf_builtin -I../lhapdf

########################################################################
# For the moment, the submodule dependencies will be hard-coded
prc_core_sub.lo: prc_core.lo
prc_test_core_sub.lo: prc_test_core.lo
prc_template_me_sub.lo: prc_template_me.lo
prc_omega_sub.lo: prc_omega.lo
prc_external_sub.lo: prc_external.lo
prc_threshold_sub.lo: prc_threshold.lo

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

PRELUDE = $(top_srcdir)/src/noweb-frame/whizard-prelude.nw
POSTLUDE = $(top_srcdir)/src/noweb-frame/whizard-postlude.nw

me_methods.stamp: $(PRELUDE) $(srcdir)/me_methods.nw $(POSTLUDE)
	@rm -f me_methods.tmp
	@touch me_methods.tmp
	for src in \
	  $(libme_methods_la_SOURCES) \
	  $(libme_methods_ut_la_SOURCES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	@mv -f me_methods.tmp me_methods.stamp

$(libme_methods_la_SOURCES) $(libme_methods_ut_la_SOURCES): me_methods.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f me_methods.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) me_methods.stamp; \
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
	-rm -f me_methods.stamp me_methods.tmp
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
