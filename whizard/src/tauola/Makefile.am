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

## The files in this directory end up in an auxiliary libtool library.
noinst_LTLIBRARIES = libtauola_interface.la
libtauola_interface_la_SOURCES = \
  $(TAUOLA_MODULES) \
  $(TAUOLA_SUBMODULES)

TAUOLA_MODULES = \
  hepev4_aux.f90 \
  tauola_interface.f90

TAUOLA_SUBMODULES = \
  tauola_interface_sub.f90

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ${TAUOLA_MODULES:.f90=.$(FCMOD)}

libtauola_interface_Modules = ${TAUOLA_MODULES:.f90=}
Modules: Makefile
	@for module in $(libtauola_interface_Modules); do \
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
  ../system/Modules \
  ../variables/Modules \
  ../qft/Modules \
  ../events/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(libtauola_interface_la_SOURCES)
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
Makefile.depend: $(libtauola_interface_la_SOURCES)
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

if PYTHIA6_AVAILABLE
libtauola_interface_la_LIBADD = ../../pythia6/libpythia6_wo.la ../../tauola/libtauola_wo.la
else
libtauola_interface_la_LIBADD = ../../pythia6/libpythia6_wo_dummy.la ../../tauola/libtauola_wo_dummy.la
endif

AM_FCFLAGS =

########################################################################
# For the moment, the submodule dependencies will be hard-coded
tauola_interface_sub.lo: tauola_interface.lo

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

########################################################################
## Explicit dependencies

AM_FCFLAGS += -I../basics -I../utilities -I../system -I../expr_base -I../combinatorics -I../parsing -I../physics -I../qft -I../types -I../fastjet -I../particles -I../events -I../variables

MODULES= \
  tauola_interface.$(FCMOD) \
  hepev4_aux.$(FCMOD)

########################################################################
## Non-standard cleanup tasks

## Remove backup files
maintainer-clean-local:
	-rm -f *~

## Remove module files
clean-local:
	-rm -f $(MODULES)
if FC_SUBMODULES
	-rm -f *.smod *.sub
endif
