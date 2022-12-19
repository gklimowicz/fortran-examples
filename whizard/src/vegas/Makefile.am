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

## The files in this directory implement the VEGAS algorithm and the
## multi-channel implementation with VEGAS as backbone integrator

## We create a library which is still to be combined with auxiliary libs.
noinst_LTLIBRARIES = libvegas.la
check_LTLIBRARIES = libvegas_ut.la

COMMON_F90 =

MPI_F90 = \
  $(MPI_MODULES) \
  $(MPI_SUBMODULES)

MPI_MODULES = \
  vegas.f90_mpi \
  vamp2.f90_mpi
MPI_SUBMODULES = \
  vegas_sub.f90_mpi \
  vamp2_sub.f90_mpi

LOAD_BALANCER_MODULES = \
  request_base.f90 \
  request_simple.f90 \
  request_caller.f90 \
  request_state.f90 \
  balancer_base.f90 \
  balancer_simple.f90 \
  balancer_channel.f90 \
  request_callback.f90

LOAD_BALANCER_SUBMODULES = \
  request_base_sub.f90 \
  request_simple_sub.f90 \
  request_caller_sub.f90 \
  request_state_sub.f90 \
  balancer_base_sub.f90 \
  balancer_simple_sub.f90 \
  balancer_channel_sub.f90 \
  request_callback_sub.f90

SERIAL_F90 = \
  $(SERIAL_MODULES) \
  $(SERIAL_SUBMODULES)

SERIAL_MODULES = \
  vegas.f90_serial \
  vamp2.f90_serial
SERIAL_SUBMODULES = \
  vegas_sub.f90_serial \
  vamp2_sub.f90_serial

EXTRA_DIST = \
  $(COMMON_F90) \
  $(SERIAL_F90) \
  $(MPI_F90) \
  $(LOAD_BALANCER_MODULES) \
  $(LOAD_BALANCER_SUBMODULES)

nodist_libvegas_la_SOURCES = \
  $(COMMON_F90) \
  $(VEGAS_MODULES) \
  $(VEGAS_SUBMODULES)

VEGAS_MODULES = \
  vegas.f90 \
  vamp2.f90

VEGAS_SUBMODULES = \
  vegas_sub.f90 \
  vamp2_sub.f90

if FC_USE_MPI
VEGAS_MODULES += \
  $(LOAD_BALANCER_MODULES)
VEGAS_SUBMODULES += \
  $(LOAD_BALANCER_SUBMODULES)
endif

DISTCLEANFILES = \
  vegas.f90 \
  vegas_sub.f90 \
  vamp2.f90 \
  vamp2_sub.f90

if FC_USE_MPI
vegas.f90: vegas.f90_mpi
	-cp -f $< $@
vegas_sub.f90: vegas_sub.f90_mpi
	-cp -f $< $@
vamp2.f90: vamp2.f90_mpi
	-cp -f $< $@
vamp2_sub.f90: vamp2_sub.f90_mpi
	-cp -f $< $@
else
vegas.f90: vegas.f90_serial
	-cp -f $< $@
vegas_sub.f90: vegas_sub.f90_serial
	-cp -f $< $@
vamp2.f90: vamp2.f90_serial
	-cp -f $< $@
vamp2_sub.f90: vamp2_sub.f90_serial
	-cp -f $< $@
endif

libvegas_ut_la_SOURCES = \
  vegas_uti.f90 vegas_ut.f90 \
  vamp2_uti.f90 vamp2_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = vegas.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ${VEGAS_MODULES:.f90=.$(FCMOD)}

libvegas_Modules = $(VEGAS_MODULES:.f90=) $(libvegas_ut_la_SOURCES:.f90=)
Modules: Makefile
	@for module in $(libvegas_Modules); do \
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
  ../rng/Modules \
  ../phase_space/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(nodist_libvegas_la_SOURCES) $(libvegas_ut_la_SOURCES)
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
		done ; \
	done

DISTCLEANFILES += Module_dependencies.sed


# The following line just says
#    include Makefile.depend
# but in a portable fashion (depending on automake's AM_MAKE_INCLUDE
@am__include@ @am__quote@Makefile.depend@am__quote@

Makefile.depend: Module_dependencies.sed
Makefile.depend: $(nodist_libvegas_la_SOURCES) $(libvegas_ut_la_SOURCES)
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

AM_FCFLAGS = -I../basics -I../utilities -I../testing -I../system -I../combinatorics -I../rng -I../physics -I../fastjet -I../qft -I../matrix_elements -I../types -I../particles -I../beams -I../phase_space -I../expr_base -I../variables -I../../vamp/src -I../mci

########################################################################
vegas_sub.lo: vegas.lo
vamp2_sub.lo: vamp2.lo
if FC_USE_MPI
request_base_sub.lo: request_base.lo
request_simple_sub.lo: request_simple.lo
request_caller_sub.lo: request_caller.lo
request_state_sub.lo: request_state.lo
balancer_base_sub.lo: balancer_base.lo
balancer_simple_sub.lo: balancer_simple.lo
balancer_channel_sub.lo: balancer_channel.lo
request_callback_sub.lo: request_callback.lo
endif

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

NOMPI_FILTER = -filter "sed 's/defn NOMPI:/defn/'"
MPI_FILTER = -filter "sed 's/defn MPI:/defn/'"

PRELUDE = $(top_srcdir)/src/noweb-frame/whizard-prelude.nw
POSTLUDE = $(top_srcdir)/src/noweb-frame/whizard-postlude.nw

vegas.stamp: $(PRELUDE) $(srcdir)/vegas.nw $(POSTLUDE)
	@rm -f vegas.tmp
	@touch vegas.tmp
	for src in $(COMMON_F90) $(libvegas_ut_la_SOURCES); do \
		$(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
		done
	for src in $(SERIAL_F90:.f90_serial=.f90); do \
		$(NOTANGLE) -R[[$$src]] $(NOMPI_FILTER) $^ | $(CPIF) $$src'_serial'; \
		done
	for src in $(MPI_F90:.f90_mpi=.f90); do \
		$(NOTANGLE) -R[[$$src]] $(MPI_FILTER) $^ | $(CPIF) $$src'_mpi'; \
		done
	for src in $(LOAD_BALANCER_MODULES) $(LOAD_BALANCER_SUBMODULES); do \
		$(NOTANGLE) -R[[$$src]] $(MPI_FILTER) $^ | $(CPIF) $$src; \
		done
	@mv -f vegas.tmp vegas.stamp

$(COMMON_F90) $(SERIAL_F90) $(MPI_F90) $(LOAD_BALANCER_MODULES) $(LOAD_BALANCER_SUBMODULES) $(libvegas_ut_la_SOURCES): vegas.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
		rm -f vegas.stamp; \
		$(MAKE) $(AM_MAKEFLAGS) vegas.stamp; \
	fi

endif


########################################################################
## Non-standard cleanup tasks
## Remove sources that can be recreated using NOWEB
if NOWEB_AVAILABLE
maintainer-clean-noweb:
	-rm -f *.f90 *.f90_serial *.f90_mpi *.c
endif
.PHONY: maintainer-clean-noweb

## Remove those sources also if builddir and srcdir are different
if NOWEB_AVAILABLE
clean-noweb:
	test "$(srcdir)" != "." && rm -f *.f90 *.f90_serial *.f90_mpi *.c || true
endif
.PHONY: clean-noweb

## Remove F90 module files
clean-local: clean-noweb
	-rm -f vegas.stamp vegas.tmp
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
