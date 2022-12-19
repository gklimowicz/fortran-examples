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

## If (parts of) LHAPDF are not available, we create a dummy to stub them.
if !LHAPDF5_AVAILABLE
if !LHAPDF6_AVAILABLE
noinst_LTLIBRARIES = libLHAPDF5_dummy.la
libLHAPDF5_dummy_la_SOURCES = \
	lhapdf5_full_dummy.f90 \
	lhapdf5_has_photon_dummy.f90
endif
else
if LHAPDF5_HAS_PHOTON_DUMMY
noinst_LTLIBRARIES = libLHAPDF5_dummy.la
libLHAPDF5_dummy_la_SOURCES = \
	lhapdf5_has_photon_dummy.f90
endif
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
## Non-standard cleanup tasks

## Remove backup files
maintainer-clean-local:
	-rm -f *~
