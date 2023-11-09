## Makefile.am -- Makefile for WHIZARD
##
## Process this file with automake to produce Makefile.in
##
########################################################################
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
noinst_LTLIBRARIES = libbasics.la

nodist_libbasics_la_SOURCES = \
  kinds.f90
libbasics_la_SOURCES = \
  iso_varying_string.f90 \
  io_units.f90 \
  constants.f90

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ${nodist_libbasics_la_SOURCES:.f90=.$(FCMOD)} \
  ${libbasics_la_SOURCES:.f90=.$(FCMOD)}

libbasics_Modules = \
  ${nodist_libbasics_la_SOURCES:.f90=} \
  ${libbasics_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(libbasics_Modules); do \
          echo $$module >> $@.new; \
        done
	@if diff $@ $@.new -q >/dev/null; then \
          rm $@.new; \
        else \
          mv $@.new $@; echo "Modules updated"; \
        fi
BUILT_SOURCES = Modules


SUFFIXES = .lo .$(FCMOD)
# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
	@:
#	touch $@

# Explicit dependencies
constants.lo: kinds.lo


AM_FCFLAGS =

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

.PHONY: clean-noweb

clean-local: clean-noweb
	-rm -f *.$(FCMOD)
if FC_SUBMODULES
	-rm -f *.smod *.sub
endif

DISTCLEANFILES = kinds.f90

## Remove backup files
maintainer-clean-local:
	-rm -f *~
