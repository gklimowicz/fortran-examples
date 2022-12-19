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

## The files in this directory make up the WHIZARD core

## We create a library which is still to be combined with auxiliary libs.
noinst_LTLIBRARIES = libwhizard_core.la
check_LTLIBRARIES = libwhizard_core_ut.la

COMMON_F90 = \
  user_files.f90 \
  rt_data.f90  \
  dispatch_me_methods.f90 \
  process_configurations.f90 \
  compilations.f90 \
  integrations.f90 \
  event_streams.f90 \
  restricted_subprocesses.f90 \
  commands.f90 \
  whizard.f90 \
  features.f90
MPI_F90 = \
  integrations_sub.f90_mpi \
  simulations.f90_mpi \
  simulations_sub.f90_mpi \
  commands_sub.f90_mpi
SERIAL_F90 = \
  integrations_sub.f90_serial \
  simulations.f90_serial \
  simulations_sub.f90_serial \
  commands_sub.f90_serial
COMMON_SUBMODULES = \
  user_files_sub.f90 \
  rt_data_sub.f90 \
  dispatch_me_methods_sub.f90 \
  process_configurations_sub.f90 \
  compilations_sub.f90 \
  event_streams_sub.f90 \
  restricted_subprocesses_sub.f90 \
  whizard_sub.f90 \
  features_sub.f90
WHIZARD_CORE_SUBMODULES = \
  $(COMMON_SUBMODULES) \
  integrations_sub.f90 \
  simulations_sub.f90 \
  commands_sub.f90
WHIZARD_CORE_MODULES = \
  $(COMMON_F90) \
  simulations.f90

EXTRA_DIST = \
  $(COMMON_F90) \
  $(COMMON_SUBMODULES) \
  $(SERIAL_F90) \
  $(MPI_F90)

nodist_libwhizard_core_la_SOURCES = \
  $(WHIZARD_CORE_MODULES) \
  $(WHIZARD_CORE_SUBMODULES)

DISTCLEANFILES = integrations_sub.f90 simulations.f90 simulations_sub.f90 commands_sub.f90

if FC_USE_MPI
integrations_sub.f90: integrations_sub.f90_mpi
	-cp -f $< $@
simulations.f90: simulations.f90_mpi
	-cp -f $< $@
simulations_sub.f90: simulations_sub.f90_mpi
	-cp -f $< $@
commands_sub.f90: commands_sub.f90_mpi
	-cp -f $< $@
else
integrations_sub.f90: integrations_sub.f90_serial
	-cp -f $< $@
simulations.f90: simulations.f90_serial
	-cp -f $< $@
simulations_sub.f90: simulations_sub.f90_serial
	-cp -f $< $@
commands_sub.f90: commands_sub.f90_serial
	-cp -f $< $@
endif

libwhizard_core_ut_la_SOURCES = \
  expr_tests_uti.f90 expr_tests_ut.f90 \
  rt_data_uti.f90 rt_data_ut.f90 \
  dispatch_uti.f90 dispatch_ut.f90 \
  process_configurations_uti.f90 process_configurations_ut.f90 \
  compilations_uti.f90 compilations_ut.f90 \
  integrations_uti.f90 integrations_ut.f90 \
  event_streams_uti.f90 event_streams_ut.f90 \
  restricted_subprocesses_uti.f90 restricted_subprocesses_ut.f90 \
  simulations_uti.f90 simulations_ut.f90 \
  commands_uti.f90 commands_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = whizard.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ${WHIZARD_CORE_MODULES:.f90=.$(FCMOD)} \
  ${libwhizard_core_la_SOURCES:.f90=.$(FCMOD)}

libwhizard_core_Modules = \
  ${WHIZARD_CORE_MODULES:.f90=} \
  ${libwhizard_core_ut_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(libwhizard_core_Modules); do \
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
  ../types/Modules \
  ../matrix_elements/Modules \
  ../particles/Modules \
  ../beams/Modules \
  ../me_methods/Modules \
  ../pythia8/Modules \
  ../events/Modules \
  ../phase_space/Modules \
  ../mci/Modules \
  ../vegas/Modules \
  ../blha/Modules \
  ../gosam/Modules \
  ../openloops/Modules \
  ../recola/Modules \
  ../fks/Modules \
  ../variables/Modules \
  ../model_features/Modules \
  ../muli/Modules \
  ../shower/Modules \
  ../matching/Modules \
  ../process_integration/Modules \
  ../transforms/Modules \
  ../threshold/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(nodist_libwhizard_core_la_SOURCES) \
  $(libwhizard_core_ut_la_SOURCES) \
  $(EXTRA_libwhizard_core_la_SOURCES)
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
Makefile.depend: $(nodist_libwhizard_core_la_SOURCES) \
  $(libwhizard_core_ut_la_SOURCES) \
  $(EXTRA_libwhizard_core_la_SOURCES)
#  $(nodist_libwhizard_main_la_SOURCES) \
#  $(nodist_libwhizard_main_ut_la_SOURCES) \
	@rm -f $@
	for src in $^; do \
		module="`basename $$src | sed 's/\.f[90][0358]//'`"; \
		grep '^ *use ' $$src \
			| grep -v '!NODEP!' \
			| sed -e 's/^ *use */'$$module'.lo: /' \
			-f Module_dependencies.sed; \
	done > $@

DISTCLEANFILES += Makefile.depend

SUFFIXES = .lo .$(FCMOD)
# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
	@:
#	touch $@

AM_FCFLAGS = -I../basics -I../utilities -I../testing -I../system -I../combinatorics -I../parsing -I../rng -I../physics -I../qed_pdf -I../qft -I../expr_base -I../types -I../matrix_elements -I../particles -I../beams -I../me_methods -I../events -I../phase_space -I../mci -I../vegas  -I../blha -I../gosam -I../openloops -I../fks -I../variables -I../model_features -I../muli -I../pythia8 -I../shower -I../matching -I../process_integration -I../transforms -I../xdr -I../../vamp/src -I../pdf_builtin -I../../circe1/src -I../../circe2/src -I../lhapdf -I../fastjet -I../threshold -I../tauola -I../recola

########################################################################
user_files_sub.lo: user_files.lo
rt_data_sub.lo: rt_data.lo
dispatch_me_methods_sub.lo: dispatch_me_methods.lo
process_configurations_sub.lo: process_configurations.lo
compilations_sub.lo: compilations.lo
event_streams_sub.lo: event_streams.lo
integrations_sub.lo: integrations.lo
restricted_subprocesses_sub.lo: restricted_subprocesses.lo
simulations_sub.lo: simulations.lo
commands_sub.lo: commands.lo
whizard_sub.lo: whizard.lo
features_sub.lo: features.lo

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

FILTER = -filter "sed 's/defn MPI:/defn/'"

COMMON_SRC = \
  $(COMMON_F90) \
  $(libwhizard_core_ut_la_SOURCES) \
  $(EXTRA_libwhizard_core_la_SOURCES)

PRELUDE = $(top_srcdir)/src/noweb-frame/whizard-prelude.nw
POSTLUDE = $(top_srcdir)/src/noweb-frame/whizard-postlude.nw

whizard.stamp: $(PRELUDE) $(srcdir)/whizard.nw $(POSTLUDE)
	@rm -f whizard.tmp
	@touch whizard.tmp
	for src in $(COMMON_SRC); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	for src in $(COMMON_SUBMODULES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
	done
	for src in $(MPI_F90:.f90_mpi=.f90); do \
	  $(NOTANGLE) -R[[$$src]] $(FILTER) $^ | $(CPIF) $$src'_mpi'; \
        done
	for src in $(SERIAL_F90:.f90_serial=.f90); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src'_serial'; \
        done
	@mv -f whizard.tmp whizard.stamp

$(MPI_F90) $(SERIAL_F90) $(COMMON_SRC) $(COMMON_SUBMODULES): whizard.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f whizard.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) whizard.stamp; \
	fi

endif

########################################################################
## Non-standard cleanup tasks
## Remove sources that can be recreated using NOWEB
if NOWEB_AVAILABLE
maintainer-clean-noweb:
	-rm -f *.f90 *.f90_mpi *.f90_serial *.c
endif
.PHONY: maintainer-clean-noweb

## Remove those sources also if builddir and srcdir are different
if NOWEB_AVAILABLE
clean-noweb:
	test "$(srcdir)" != "." && rm -f *.f90 *.f90_mpi *.f90_serial *.c || true
endif
.PHONY: clean-noweb

## Remove F90 module files
clean-local: clean-noweb
	-rm -f whizard.stamp whizard.tmp
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
