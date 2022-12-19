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
noinst_LTLIBRARIES = libapi.la
check_LTLIBRARIES = libapi_ut.la
check_LTLIBRARIES += libapi_ut_c.la
check_LTLIBRARIES += libapi_ut_cc.la

COMMON_F90 = \
  api.f90 \
  api_c.f90

COMMON_C =
MPI_C =
SERIAL_C =

COMMON_CC = \
  api_cc.cc
MPI_CC =
SERIAL_CC =

EXTRA_DIST = \
  $(COMMON_C) \
  $(SERIAL_C) \
  $(MPI_C) \
  $(COMMON_CC) \
  $(SERIAL_CC) \
  $(MPI_CC)

libapi_la_SOURCES = \
  $(COMMON_F90) \
  api_cc.cc

DISTCLEANFILES = api.f90

libapi_ut_la_SOURCES = \
  api_uti.f90 api_ut.f90 \
  api_hepmc_uti.f90 api_hepmc_ut.f90 \
  api_lcio_uti.f90 api_lcio_ut.f90


nodist_libapi_ut_c_la_SOURCES = \
  api_ut_c.c

DISTCLEANFILES += api_ut_c.c

MPI_C += api_ut_c.c_mpi
SERIAL_C += api_ut_c.c_serial

if FC_USE_MPI
api_ut_c.c: api_ut_c.c_mpi
	-cp -f $< $@
else
api_ut_c.c: api_ut_c.c_serial
	-cp -f $< $@
endif

nodist_libapi_ut_cc_la_SOURCES = \
  whizard_ut.cc \
  api_ut_cc.cc

DISTCLEANFILES += \
  whizard_ut.cc \
  api_ut_cc.cc

COMMON_CC += whizard_ut.cc
MPI_CC += api_ut_cc.cc_mpi
SERIAL_CC += api_ut_cc.cc_serial

if FC_USE_MPI
api_ut_cc.cc: api_ut_cc.cc_mpi
	-cp -f $< $@
else
api_ut_cc.cc: api_ut_cc.cc_serial
	-cp -f $< $@
endif



libapi_la_CPPFLAGS =
libapi_ut_cc_la_CPPFLAGS =
if HEPMC3_AVAILABLE
libapi_la_CPPFLAGS += -DWHIZARD_WITH_HEPMC3 $(HEPMC_INCLUDES)
libapi_ut_cc_la_CPPFLAGS += -DWHIZARD_WITH_HEPMC3 $(HEPMC_INCLUDES)
endif
if HEPMC2_AVAILABLE
libapi_la_CPPFLAGS += -DWHIZARD_WITH_HEPMC2 $(HEPMC_INCLUDES)
libapi_ut_cc_la_CPPFLAGS += -DWHIZARD_WITH_HEPMC2 $(HEPMC_INCLUDES)
endif
if LCIO_AVAILABLE
libapi_la_CPPFLAGS += -DWHIZARD_WITH_LCIO $(LCIO_INCLUDES)
libapi_ut_cc_la_CPPFLAGS += -DWHIZARD_WITH_LCIO $(LCIO_INCLUDES)
endif

include_HEADERS = \
  whizard.h
dist_noinst_HEADERS = \
  whizard_ut.h

## Omitting this would exclude it from the distribution
dist_noinst_DATA = api.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  api.$(FCMOD)

# Dump module names into file Modules
libapi_Modules = \
  ${libapi_la_SOURCES:.f90=} \
  ${libapi_ut_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(libapi_Modules); do \
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
  ../threshold/Modules \
  ../whizard-core/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(libapi_la_SOURCES) \
  $(libapi_ut_la_SOURCES)
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
  $(libapi_la_SOURCES) \
  $(libapi_ut_la_SOURCES)
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

AM_FCFLAGS = -I../basics -I../utilities -I../testing -I../system -I../combinatorics -I../parsing -I../rng -I../physics -I../qed_pdf -I../qft -I../expr_base -I../types -I../matrix_elements -I../particles -I../beams -I../me_methods -I../events -I../phase_space -I../mci -I../vegas -I../blha -I../gosam -I../openloops -I../fks -I../variables -I../model_features -I../muli -I../pythia8 -I../shower -I../matching -I../process_integration -I../transforms -I../xdr -I../../vamp/src -I../pdf_builtin -I../../circe1/src -I../../circe2/src -I../lhapdf -I../fastjet -I../threshold -I../tauola -I../recola -I../whizard-core
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
  $(COMMON_CC) \
  $(libapi_ut_la_SOURCES) \
  $(libapi_ut_cc_la_SOURCES) \
  $(include_HEADERS) \
  $(dist_noinst_HEADERS)

PRELUDE = $(top_srcdir)/src/noweb-frame/whizard-prelude.nw
POSTLUDE = $(top_srcdir)/src/noweb-frame/whizard-postlude.nw

api.stamp: $(PRELUDE) $(srcdir)/api.nw $(POSTLUDE)
	@rm -f api.tmp
	@touch api.tmp
	for src in $(COMMON_SRC); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	for src in $(MPI_C:.c_mpi=.c); do \
	  $(NOTANGLE) -R[[$$src]] $(FILTER) $^ | $(CPIF) $$src'_mpi'; \
        done
	for src in $(SERIAL_C:.c_serial=.c); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src'_serial'; \
        done
	for src in $(MPI_CC:.cc_mpi=.cc); do \
	  $(NOTANGLE) -R[[$$src]] $(FILTER) $^ | $(CPIF) $$src'_mpi'; \
        done
	for src in $(SERIAL_CC:.cc_serial=.cc); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src'_serial'; \
        done
	@mv -f api.tmp api.stamp

$(COMMON_SRC) $(MPI_C) $(SERIAL_C) $(MPI_CC) $(SERIAL_CC): api.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f api.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) api.stamp; \
	fi

endif

########################################################################
## Non-standard cleanup tasks
## Remove sources that can be recreated using NOWEB
if NOWEB_AVAILABLE
maintainer-clean-noweb:
	-rm -f *.f90 *.f90_mpi *.f90_serial *.c *.cc *.h *.cpp
endif
.PHONY: maintainer-clean-noweb

## Remove those sources also if builddir and srcdir are different
if NOWEB_AVAILABLE
clean-noweb:
	test "$(srcdir)" != "." && rm -f *.f90 *.f90_mpi *.f90_serial \
	   *.c_mpi *.c_serial *.cc_mpi *.cc_serial *.c *.cc *.h *.cpp || true
endif
.PHONY: clean-noweb

## Remove F90 module files
clean-local: clean-noweb
	-rm -f api.stamp api.tmp
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
