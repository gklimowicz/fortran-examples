## Makefile.am -- Makefile for WHIZARD
##
## Process this file with automake to produce Makefile.in
#
# Copyright (C) 1999-2022 by
#     Wolfgang Kilian <kilian@physik.uni-siegen.de>
#     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
#     Juergen Reuter <juergen.reuter@desy.de>
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

## The files in this directory form the core of process integration

## We create a library which is still to be combined with auxiliary libs.
noinst_LTLIBRARIES = libprocess_integration.la
check_LTLIBRARIES = libprocess_integration_ut.la

libprocess_integration_la_SOURCES = \
  $(PROCINT_MODULES) \
  $(PROCINT_SUBMODULES)

PROCINT_MODULES = \
  subevt_expr.f90 \
  parton_states.f90 \
  pcm_base.f90 \
  pcm.f90 \
  process_counter.f90 \
  process_config.f90 \
  process_mci.f90 \
  process.f90 \
  kinematics.f90 \
  instances.f90 \
  process_stacks.f90

PROCINT_SUBMODULES = \
  subevt_expr_sub.f90 \
  parton_states_sub.f90 \
  pcm_base_sub.f90 \
  pcm_sub.f90 \
  process_counter_sub.f90 \
  process_config_sub.f90 \
  process_mci_sub.f90 \
  process_sub.f90 \
  kinematics_sub.f90 \
  instances_sub.f90 \
  process_stacks_sub.f90

libprocess_integration_ut_la_SOURCES = \
  processes_uti.f90 processes_ut.f90 \
  process_stacks_uti.f90 process_stacks_ut.f90 \
  parton_states_uti.f90 parton_states_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = process_integration.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ${PROCINT_MODULES:.f90=.$(FCMOD)}

libprocess_integration_Modules = \
  ${PROCINT_MODULES:.f90=} \
  ${libprocess_integration_ut_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(libprocess_integration_Modules); do \
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
  ../physics/Modules \
  ../qft/Modules \
  ../expr_base/Modules \
  ../types/Modules \
  ../matrix_elements/Modules \
  ../particles/Modules \
  ../beams/Modules \
  ../me_methods/Modules \
  ../model_features/Modules \
  ../phase_space/Modules \
  ../mci/Modules \
  ../blha/Modules \
  ../gosam/Modules \
  ../openloops/Modules \
  ../recola/Modules \
  ../fks/Modules \
  ../variables/Modules \
  ../parsing/Modules \
  ../shower/Modules \
  ../threshold/Modules

include_modules_bare = ${module_lists:/Modules=}
include_modules = ${include_modules_bare:../%=-I../%}

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(libprocess_integration_la_SOURCES) \
  $(libprocess_integration_ut_la_SOURCES)
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
Makefile.depend: $(libprocess_integration_la_SOURCES) \
  $(libprocess_integration_ut_la_SOURCES)
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

AM_FCFLAGS = $(include_modules) -I../../vamp/src -I../../circe1/src -I../../circe2/src -I../pdf_builtin -I../lhapdf -I../qed_pdf -I../fastjet

########################################################################
subevt_expr_sub.lo: subevt_expr.lo
parton_states_sub.lo: parton_states.lo
pcm_base_sub.lo: pcm_base.lo
pcm_sub.lo: pcm.lo
process_counter_sub.lo: process_counter.lo
process_config_sub.lo: process_config.lo
process_mci_sub.lo: process_mci.lo
process_sub.lo: process.lo
kinematics_sub.lo: kinematics.lo
instances_sub.lo: instances.lo
process_stacks_sub.lo: process_stacks.lo

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

if RECOLA_AVAILABLE
AM_FCFLAGS += $(RECOLA_INCLUDES)
endif

########################################################################
## Non-standard targets and dependencies

## (Re)create F90 sources from NOWEB source.
if NOWEB_AVAILABLE

PRELUDE = $(top_srcdir)/src/noweb-frame/whizard-prelude.nw
POSTLUDE = $(top_srcdir)/src/noweb-frame/whizard-postlude.nw

process_integration.stamp: $(PRELUDE) $(srcdir)/process_integration.nw $(POSTLUDE)
	@rm -f process_integration.tmp
	@touch process_integration.tmp
	for src in \
          $(libprocess_integration_la_SOURCES) \
          $(libprocess_integration_ut_la_SOURCES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	@mv -f process_integration.tmp process_integration.stamp

$(libprocess_integration_la_SOURCES) $(libprocess_integration_ut_la_SOURCES): process_integration.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f process_integration.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) process_integration.stamp; \
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
	-rm -f process_integration.stamp process_integration.tmp
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
