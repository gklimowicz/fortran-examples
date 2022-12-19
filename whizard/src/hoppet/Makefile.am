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
if HOPPET_AVAILABLE
noinst_LTLIBRARIES = libhoppet.la
libhoppet_la_SOURCES = hoppet.f90
libhoppet_la_FCFLAGS = $(AM_FCFLAGS) $(HOPPET_INCLUDES) 
else
noinst_LTLIBRARIES = libhoppet_dummy.la
libhoppet_dummy_la_SOURCES = hoppet_dummy.f90
endif

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

libhoppet.lo: \
  ../system/system_dependencies.$(FCMOD) \
  ../pdf_builtin/pdf_builtin.$(FCMOD) \
  ../lhapdf/lhapdf.$(FCMOD)

libhoppet_dummy.lo: \
  ../lhapdf/lhapdf.$(FCMOD)

AM_FCFLAGS += -I../basics -I../system -I../pdf_builtin -I../lhapdf

MODULES= \
   hoppet.$(FCMOD)

########################################################################
## Cleanup

## Remove backup files
maintainer-clean-local:
	-rm -f *~
