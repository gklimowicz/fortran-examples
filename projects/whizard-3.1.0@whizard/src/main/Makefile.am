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

## The main program is contained in a library on its own.
lib_LTLIBRARIES = libwhizard_main.la
check_LTLIBRARIES = libwhizard_main_ut.la

COMMON_F90 = \
  cmdline_options.f90
MPI_F90 = \
  main.f90_mpi
SERIAL_F90 = \
  main.f90_serial

nodist_libwhizard_main_la_SOURCES = \
  $(COMMON_F90) \
  main.f90

DISTCLEANFILES = main.f90

if FC_USE_MPI
main.f90: main.f90_mpi
	-cp -f $< $@
else
main.f90: main.f90_serial
	-cp -f $< $@
endif

MPI_F90 += main_ut.f90_mpi
SERIAL_F90 += main_ut.f90_serial

nodist_libwhizard_main_ut_la_SOURCES = \
  $(COMMON_F90) \
  main_ut.f90

DISTCLEANFILES += main_ut.f90

if FC_USE_MPI
main_ut.f90: main_ut.f90_mpi
	-cp -f $< $@
else
main_ut.f90: main_ut.f90_serial
	-cp -f $< $@
endif

EXTRA_DIST = \
  $(COMMON_F90) \
  $(SERIAL_F90) \
  $(MPI_F90)

## Omitting this would exclude it from the distribution
dist_noinst_DATA = main.nw

# # Dump module names into file Modules
# libwhizard_core_Modules = \
#   ${libwhizard_core_la_SOURCES:.f90=} \
#   ${nodist_libwhizard_core_la_SOURCES:.f90=} \
#   ${libwhizard_core_ut_la_SOURCES:.f90=}
# Modules: Makefile
# 	@for module in $(libwhizard_core_Modules); do \
#           echo $$module >> $@.new; \
#         done
# 	@if diff $@ $@.new -q >/dev/null; then \
#           rm $@.new; \
# 	else \
#           mv $@.new $@; echo "Modules updated"; \
#         fi
# BUILT_SOURCES = Modules

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
  ../qed_pdf/Modules \
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
  ../threshold/Modules \
  ../whizard-core/Modules \
  ../api/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: \
  $(nodist_libwhizard_main_la_SOURCES) \
  $(nodist_libwhizard_main_ut_la_SOURCES)
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
Makefile.depend: \
  $(nodist_libwhizard_main_la_SOURCES) \
  $(nodist_libwhizard_main_ut_la_SOURCES)
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

AM_FCFLAGS = -I../basics -I../utilities -I../testing -I../system -I../combinatorics -I../parsing -I../rng -I../physics -I../qed_pdf -I../qft -I../expr_base -I../types -I../matrix_elements -I../particles -I../beams -I../me_methods -I../events -I../phase_space -I../mci -I../vegas -I../blha -I../gosam -I../openloops -I../fks -I../variables -I../model_features -I../muli -I../pythia8 -I../shower -I../matching -I../process_integration -I../transforms -I../xdr -I../../vamp/src -I../pdf_builtin -I../../circe1/src -I../../circe2/src -I../lhapdf -I../fastjet -I../threshold -I../tauola -I../recola -I../whizard-core -I../api
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
  $(COMMON_F90)

PRELUDE = $(top_srcdir)/src/noweb-frame/whizard-prelude.nw
POSTLUDE = $(top_srcdir)/src/noweb-frame/whizard-postlude.nw

main.stamp: $(PRELUDE) $(srcdir)/main.nw $(POSTLUDE)
	@rm -f main.tmp
	@touch main.tmp
	for src in $(COMMON_SRC); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	for src in $(MPI_F90:.f90_mpi=.f90); do \
	  $(NOTANGLE) -R[[$$src]] $(FILTER) $^ | $(CPIF) $$src'_mpi'; \
        done
	for src in $(SERIAL_F90:.f90_serial=.f90); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src'_serial'; \
        done
	@mv -f main.tmp main.stamp

$(MPI_F90) $(SERIAL_F90) $(COMMON_SRC): main.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f main.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) main.stamp; \
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
	-rm -f main.stamp main.tmp
	-rm -f *.$(FCMOD)
if FC_SUBMODULES
	-rm -f *.smod
endif

## Remove backup files
maintainer-clean-backup:
	-rm -f *~
.PHONY: maintainer-clean-backup

## Register additional clean targets
maintainer-clean-local: maintainer-clean-noweb maintainer-clean-backup
