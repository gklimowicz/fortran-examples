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
#AM_FFLAGS = -fno-automatic -ffixed-line-length-132 -fno-backslash -fPIC
AM_FFLAGS =
AM_FCFLAGS =

if PYTHIA6_AVAILABLE

if FC_IS_NAG
AM_FFLAGS += -dcfuns -w
AM_FCFLAGS += -dcfuns -w
endif

noinst_LTLIBRARIES = libtauola_wo.la

libtauola_wo_la_SOURCES = \
  tauface-jetset.f formf.f photos.f tauola.f tauola_photos_ini.f

else

noinst_LTLIBRARIES = libtauola_wo_dummy.la
libtauola_wo_dummy_la_SOURCES = tauola_dummy.f90

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

## Remove F90 module files
clean-local:
	-rm -f *.$(FC_MODULE_EXT)
if FC_SUBMODULES
	-rm -f *.smod
endif

## Remove backup files
maintainer-clean-local:
	-rm -f *~
