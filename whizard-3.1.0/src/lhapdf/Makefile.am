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
if LHAPDF6_AVAILABLE
noinst_LTLIBRARIES = libLHAPDFWrap.la
libLHAPDFWrap_la_SOURCES = \
   LHAPDFWrap.cpp lhapdf.f90
libLHAPDFWrap_la_CPPFLAGS = $(LHAPDF_CXXFLAGS)
else
noinst_LTLIBRARIES = libLHAPDFWrap_dummy.la
libLHAPDFWrap_dummy_la_SOURCES = \
   LHAPDFWrap_dummy.f90 lhapdf.f90
endif

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  lhapdf.$(FCMOD)

SUFFIXES = .lo .$(FCMOD)

# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
        @:
#       touch $@

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
## Explicit dependencies

lhapdf.lo: ../basics/kinds.$(FCMOD)
lhapdf.$(FCMOD): lhapdf.lo

AM_FCFLAGS += -I../basics

MODULES= \
  lhapdf.$(FCMOD)

########################################################################
## Non-standard cleanup tasks

## Remove backup files
maintainer-clean-local:
	-rm -f *~

## Remove module files
clean-local:
	-rm -f $(MODULES)
if FC_SUBMODULES
	-rm -f *.smod
endif
