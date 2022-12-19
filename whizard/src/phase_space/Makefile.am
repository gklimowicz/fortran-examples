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
noinst_LTLIBRARIES = libphase_space.la
check_LTLIBRARIES = libphase_space_ut.la

libphase_space_la_SOURCES = \
  $(PHS_MODULES) \
  $(PHS_SUBMODULES)

PHS_MODULES = \
  phs_base.f90 \
  phs_none.f90 \
  phs_single.f90 \
  phs_rambo.f90 \
  resonances.f90 \
  mappings.f90 \
  phs_trees.f90 \
  phs_forests.f90 \
  cascades.f90 \
  cascades2_lexer.f90 \
  cascades2.f90 \
  phs_wood.f90 \
  phs_fks.f90 \
  dispatch_phase_space.f90

PHS_SUBMODULES = \
  phs_base_sub.f90 \
  phs_none_sub.f90 \
  phs_single_sub.f90 \
  phs_rambo_sub.f90 \
  resonances_sub.f90 \
  mappings_sub.f90 \
  phs_trees_sub.f90 \
  phs_forests_sub.f90 \
  cascades_sub.f90 \
  cascades2_lexer_sub.f90 \
  cascades2_sub.f90 \
  phs_wood_sub.f90 \
  phs_fks_sub.f90 \
  dispatch_phase_space_sub.f90

libphase_space_ut_la_SOURCES = \
  phs_base_uti.f90 phs_base_ut.f90 \
  phs_none_uti.f90 phs_none_ut.f90 \
  phs_single_uti.f90 phs_single_ut.f90 \
  phs_rambo_uti.f90 phs_rambo_ut.f90 \
  resonances_uti.f90 resonances_ut.f90 \
  phs_trees_uti.f90 phs_trees_ut.f90 \
  phs_forests_uti.f90 phs_forests_ut.f90 \
  cascades_uti.f90 cascades_ut.f90 \
  cascades2_lexer_uti.f90 cascades2_lexer_ut.f90 \
  cascades2_uti.f90 cascades2_ut.f90 \
  phs_wood_uti.f90 phs_wood_ut.f90 \
  phs_fks_uti.f90 phs_fks_ut.f90 \
  dispatch_phs_uti.f90 dispatch_phs_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = phase_space.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ${PHS_MODULES:.f90=.$(FCMOD)}

# Dump module names into file Modules
libphase_space_Modules = \
  ${PHS_MODULES:.f90=} \
  ${libphase_space_ut_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(libphase_space_Modules); do \
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
  ../physics/Modules \
  ../qft/Modules \
  ../types/Modules \
  ../matrix_elements/Modules \
  ../beams/Modules \
  ../model_features/Modules \
  ../variables/Modules \
  ../expr_base/Modules \
  ../threshold/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(libphase_space_la_SOURCES) $(libphase_space_ut_la_SOURCES)
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
Makefile.depend: $(libphase_space_la_SOURCES) $(libphase_space_ut_la_SOURCES)
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

AM_FCFLAGS = -I../basics -I../utilities -I../testing -I../system -I../combinatorics -I../parsing -I../physics -I../fastjet -I../qed_pdf -I../qft -I../matrix_elements -I../types -I../particles -I../beams -I../rng -I../../circe1/src -I../../circe2/src  -I../pdf_builtin -I../lhapdf -I../model_features -I../variables -I../expr_base -I../threshold

########################################################################
phs_base_sub.lo: phs_base.lo
phs_none_sub.lo: phs_none.lo
phs_single_sub.lo: phs_single.lo
phs_rambo_sub.lo: phs_rambo.lo
resonances_sub.lo: resonances.lo
mappings_sub.lo: mappings.lo
phs_trees_sub.lo: phs_trees.lo
phs_forests_sub.lo: phs_forests.lo
cascades_sub.lo: cascades.lo
cascades2_lexer_sub.lo: cascades2_lexer.lo
cascades2_sub.lo: cascades2.lo
phs_wood_sub.lo: phs_wood.lo
phs_fks_sub.lo: phs_fks.lo
dispatch_phase_space_sub.lo: dispatch_phase_space.lo

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

phase_space.stamp: $(PRELUDE) $(srcdir)/phase_space.nw $(POSTLUDE)
	@rm -f phase_space.tmp
	@touch phase_space.tmp
	for src in \
	  $(libphase_space_la_SOURCES) \
	  $(libphase_space_ut_la_SOURCES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	@mv -f phase_space.tmp phase_space.stamp

$(libphase_space_la_SOURCES) $(libphase_space_ut_la_SOURCES): phase_space.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f phase_space.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) phase_space.stamp; \
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
	-rm -f phase_space.stamp phase_space.tmp
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
