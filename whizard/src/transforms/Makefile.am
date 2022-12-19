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
noinst_LTLIBRARIES = libtransforms.la
check_LTLIBRARIES = libtransforms_ut.la

libtransforms_la_SOURCES = \
  $(TRANSFORMS_MODULES) \
  $(TRANSFORMS_SUBMODULES)

TRANSFORMS_MODULES = \
  event_transforms.f90 \
  hadrons.f90 \
  resonance_insertion.f90 \
  recoil_kinematics.f90 \
  isr_epa_handler.f90 \
  decays.f90 \
  tau_decays.f90 \
  shower.f90 \
  evt_nlo.f90 \
  events.f90 \
  eio_raw.f90 \
  dispatch_transforms.f90

TRANSFORMS_SUBMODULES = \
  event_transforms_sub.f90 \
  hadrons_sub.f90 \
  resonance_insertion_sub.f90 \
  recoil_kinematics_sub.f90 \
  isr_epa_handler_sub.f90 \
  decays_sub.f90 \
  tau_decays_sub.f90 \
  shower_sub.f90 \
  evt_nlo_sub.f90 \
  events_sub.f90 \
  eio_raw_sub.f90

libtransforms_ut_la_SOURCES = \
  event_transforms_uti.f90 event_transforms_ut.f90 \
  resonance_insertion_uti.f90 resonance_insertion_ut.f90 \
  recoil_kinematics_uti.f90 recoil_kinematics_ut.f90 \
  isr_epa_handler_uti.f90 isr_epa_handler_ut.f90 \
  decays_uti.f90 decays_ut.f90 \
  shower_uti.f90 shower_ut.f90 \
  events_uti.f90 events_ut.f90 \
  eio_raw_uti.f90 eio_raw_ut.f90 \
  dispatch_transforms_uti.f90 dispatch_transforms_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = transforms.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ${TRANSFORMS_MODULES:.f90=.$(FCMOD)}

libtransforms_Modules = \
  ${TRANSFORMS_MODULES:.f90=} \
  ${libtransforms_ut_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(libtransforms_Modules); do \
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
  ../physics/Modules \
  ../qft/Modules \
  ../particles/Modules \
  ../types/Modules \
  ../matrix_elements/Modules \
  ../me_methods/Modules \
  ../rng/Modules \
  ../mci/Modules \
  ../phase_space/Modules \
  ../process_integration/Modules \
  ../matching/Modules \
  ../model_features/Modules \
  ../events/Modules \
  ../shower/Modules \
  ../pythia8/Modules \
  ../expr_base/Modules \
  ../variables/Modules \
  ../beams/Modules \
  ../threshold/Modules \
  ../fks/Modules

# Explicit dependencies, not automatically generated
event_transforms.lo:  ../../vamp/src/vamp.$(FCMOD)
shower.lo:  ../lhapdf/lhapdf.$(FCMOD)
shower.lo:  ../pdf_builtin/pdf_builtin.$(FCMOD)

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(libtransforms_la_SOURCES) $(libtransforms_ut_la_SOURCES)
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
Makefile.depend: $(libtransforms_la_SOURCES) $(libtransforms_ut_la_SOURCES)
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

AM_FCFLAGS = -I../basics -I../utilities -I../testing -I../system  -I../combinatorics -I../parsing -I../physics -I../qft -I../particles -I../types -I../matrix_elements -I../me_methods -I../rng -I../mci -I../phase_space -I../process_integration -I../matching -I../model_features -I../events -I../tauola -I../shower -I../expr_base -I../variables -I../muli -I../beams -I../blha -I../threshold -I../fks -I../openloops -I../gosam -I../recola -I../../circe1/src -I../../circe2/src -I../pdf_builtin -I../lhapdf -I../qed_pdf -I../fastjet -I../../vamp/src -I../pythia8

########################################################################
event_transforms_sub.lo: event_transforms.lo
hadrons_sub.lo: hadrons.lo
resonance_insertion_sub.lo: resonance_insertion.lo
recoil_kinematics_sub.lo: recoil_kinematics.lo
isr_epa_handler_sub.lo: isr_epa_handler.lo
decays_sub.lo: decays.lo
tau_decays_sub.lo: tau_decays.lo
shower_sub.lo: shower.lo
evt_nlo_sub.lo: evt_nlo.lo
events_sub.lo: events.lo
eio_raw_sub.lo: eio_raw.lo

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

if RECOLA_AVAILABLE
AM_FCFLAGS += $(RECOLA_INCLUDES)
endif

########################################################################
## Non-standard targets and dependencies

## (Re)create F90 sources from NOWEB source.
if NOWEB_AVAILABLE

PRELUDE = $(top_srcdir)/src/noweb-frame/whizard-prelude.nw
POSTLUDE = $(top_srcdir)/src/noweb-frame/whizard-postlude.nw

transforms.stamp: $(PRELUDE) $(srcdir)/transforms.nw $(POSTLUDE)
	@rm -f transforms.tmp
	@touch transforms.tmp
	for src in \
          $(libtransforms_la_SOURCES) \
          $(libtransforms_ut_la_SOURCES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	@mv -f transforms.tmp transforms.stamp

$(libtransforms_la_SOURCES) $(libtransforms_ut_la_SOURCES): transforms.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f transforms.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) transforms.stamp; \
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
	-rm -f transforms.stamp transforms.tmp
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
