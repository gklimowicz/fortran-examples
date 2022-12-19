## Makefile.am -- Makefile for WHIZARD
# $Id: Makefile.am 1564 2010-01-21 18:19:23Z ohl $
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
#     Fabian Bach <fabian.bach@t-online.de>
#     Bijan Chokoufe <bijan.chokoufe@desy.de>    
#     Christian Speckner <cnspeckn@googlemail.com>
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
AM_FFLAGS =
AM_FCFLAGS =

if PYTHIA6_AVAILABLE

if FC_IS_NAG
AM_FFLAGS += -dcfuns -w
AM_FCFLAGS += -dcfuns -w
endif

if PYTHIA6_IS_EH
AM_FFLAGS += -DPYTHIA6_EH
endif

noinst_LTLIBRARIES = libpythia6_wo.la

if LHAPDF5_AVAILABLE
libpythia6_wo_la_SOURCES = pythia.F
else
if LHAPDF6_AVAILABLE
libpythia6_wo_la_SOURCES = pythia.F
else
libpythia6_wo_la_SOURCES = pythia.F pythia_pdf.f
endif
endif	

else 

noinst_LTLIBRARIES = libpythia6_wo_dummy.la
libpythia6_wo_dummy_la_SOURCES = pythia6_dummy.f90

endif
########################################################################
## Default Fortran compiler options

## Profiling
if FC_USE_PROFILING
AM_FFLAGS += $(FCFLAGS_PROFILING)
AM_FCFLAGS += $(FCFLAGS_PROFILING)
endif

## OpenMP
if FC_USE_OPENMP
AM_FFLAGS += $(FCFLAGS_OPENMP)
AM_FCFLAGS += $(FCFLAGS_OPENMP)
endif

########################################################################
## Non-standard cleanup tasks

## Remove backup files
maintainer-clean-local:
	-rm -f *~
