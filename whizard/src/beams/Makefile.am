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

## We create a library which is still to be combined with auxiliary libs.
noinst_LTLIBRARIES = libbeams.la
check_LTLIBRARIES = libbeams_ut.la

libbeams_la_SOURCES = \
  $(BEAMS_MODULES) \
  $(BEAMS_SUBMODULES)

BEAMS_MODULES = \
  beam_structures.f90 \
  beams.f90 \
  sf_mappings.f90 \
  sf_aux.f90 \
  sf_base.f90 \
  sf_isr.f90 \
  sf_epa.f90 \
  sf_ewa.f90 \
  sf_escan.f90 \
  sf_gaussian.f90 \
  sf_beam_events.f90 \
  sf_circe1.f90 \
  sf_circe2.f90 \
  hoppet_interface.f90 \
  sf_pdf_builtin.f90 \
  sf_lhapdf.f90 \
  pdf.f90 \
  dispatch_beams.f90

BEAMS_SUBMODULES = \
  beam_structures_sub.f90 \
  beams_sub.f90 \
  sf_mappings_sub.f90 \
  sf_aux_sub.f90 \
  sf_base_sub.f90 \
  sf_isr_sub.f90 \
  sf_epa_sub.f90 \
  sf_ewa_sub.f90 \
  sf_escan_sub.f90 \
  sf_gaussian_sub.f90 \
  sf_beam_events_sub.f90 \
  sf_circe1_sub.f90 \
  sf_circe2_sub.f90 \
  sf_pdf_builtin_sub.f90 \
  sf_lhapdf_sub.f90 \
  pdf_sub.f90 \
  dispatch_beams_sub.f90

libbeams_ut_la_SOURCES = \
  beam_structures_uti.f90 beam_structures_ut.f90 \
  beams_uti.f90 beams_ut.f90 \
  sf_aux_uti.f90 sf_aux_ut.f90 \
  sf_mappings_uti.f90 sf_mappings_ut.f90 \
  sf_base_uti.f90 sf_base_ut.f90 \
  sf_isr_uti.f90 sf_isr_ut.f90 \
  sf_epa_uti.f90 sf_epa_ut.f90 \
  sf_ewa_uti.f90 sf_ewa_ut.f90 \
  sf_escan_uti.f90 sf_escan_ut.f90 \
  sf_gaussian_uti.f90 sf_gaussian_ut.f90 \
  sf_beam_events_uti.f90 sf_beam_events_ut.f90 \
  sf_circe1_uti.f90 sf_circe1_ut.f90 \
  sf_circe2_uti.f90 sf_circe2_ut.f90 \
  sf_pdf_builtin_uti.f90 sf_pdf_builtin_ut.f90 \
  sf_lhapdf_uti.f90 sf_lhapdf_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = beams.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ${BEAMS_MODULES:.f90=.$(FCMOD)}

libbeams_Modules = ${BEAMS_MODULES:.f90=} ${libbeams_ut_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(libbeams_Modules); do \
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
  ../types/Modules \
  ../particles/Modules \
  ../variables/Modules \
  ../qed_pdf/Modules

# Explicit dependencies, not automatically generated
sf_circe1.lo: ../../circe1/src/circe1.$(FCMOD)
sf_circe2.lo: ../../circe2/src/circe2.$(FCMOD)
hoppet_interface.lo:  ../lhapdf/lhapdf.$(FCMOD)
sf_pdf_builtin.lo: ../pdf_builtin/pdf_builtin.$(FCMOD)
sf_lhapdf.lo:  ../lhapdf/lhapdf.$(FCMOD)


$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(libbeams_la_SOURCES) $(libbeams_ut_la_SOURCES)
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
Makefile.depend: $(libbeams_la_SOURCES) $(libbeams_ut_la_SOURCES)
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

AM_FCFLAGS = -I../basics -I../utilities -I../testing -I../system -I../combinatorics -I../rng -I../physics -I../fastjet -I../qft -I../types -I../particles -I../../circe1/src -I../../circe2/src -I../pdf_builtin -I../lhapdf -I../qed_pdf -I../variables -I../expr_base -I../parsing

########################################################################
# For the moment, the submodule dependencies will be hard-coded
beam_structures_sub.lo: beam_structures.lo
beams_sub.lo: beams.lo
sf_mappings_sub.lo: sf_mappings.lo
sf_aux_sub.lo: sf_aux.lo
sf_base_sub.lo: sf_base.lo
sf_isr_sub.lo: sf_isr.lo
sf_epa_sub.lo: sf_epa.lo
sf_ewa_sub.lo: sf_ewa.lo
sf_escan_sub.lo: sf_escan.lo
sf_gaussian_sub.lo: sf_gaussian.lo
sf_beam_events_sub.lo: sf_beam_events.lo
sf_circe1_sub.lo: sf_circe1.lo
sf_circe2_sub.lo: sf_circe2.lo
sf_pdf_builtin_sub.lo: sf_pdf_builtin.lo
sf_lhapdf_sub.lo: sf_lhapdf.lo
pdf_sub.lo: pdf.lo
dispatch_beams_sub.lo: dispatch_beams.lo

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

beams.stamp: $(PRELUDE) $(srcdir)/beams.nw $(POSTLUDE)
	@rm -f beams.tmp
	@touch beams.tmp
	for src in $(libbeams_la_SOURCES) $(libbeams_ut_la_SOURCES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	@mv -f beams.tmp beams.stamp

$(libbeams_la_SOURCES) $(libbeams_ut_la_SOURCES): beams.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f beams.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) beams.stamp; \
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
	-rm -f beams.stamp beams.tmp
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
