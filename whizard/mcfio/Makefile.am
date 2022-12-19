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
AM_FCFLAGS =
AM_FFLAGS =
AM_CFLAGS = $(RPC_CFLAGS)

noinst_LTLIBRARIES = libwo_mcfio.la

libwo_mcfio_la_SOURCES = \
   mcfio_FPrintDictionary.f \
   mcfio_Util1.h \
   mcfio_Util1.c \
   mcf_ntuBldDbinc.h \
   mcf_ntuBldDbinc.c \
   mcf_NTuIOUtils.h \
   mcf_NTuIOUtils.c \
   mcf_NTuIOFiles.h \
   mcf_NTuIOFiles.c \
   mcfio_UserDictionary.h \
   mcfio_UserDictionary.c \
   mcf_evt_xdr.c \
   mcfio_FBinding.c \
   mcfio_Direct.h \
   mcfio_Direct.c \
   mcfio_SeqDummy.c \
   mcfio_Block.c \
   mcf_ntubldInit.c \
   mcf_nTupleDescript.h \
   mcf_ntubld_db.h \
   mcf_xdr.h \
   mcf_xdr_Ntuple.h \
   mcfio_Block.h \
   mcfio_Dict.h \
   mcfio_Sequential.h

SUFFIXES: .lo .$(FC_MODULE_EXT)

########################################################################
## Default Fortran compiler options

## Profiling
if FC_USE_PROFILING
AM_FCFLAGS += $(FCFLAGS_PROFILING)
AM_FFLAGS += $(FCFLAGS_PROFILING)
endif

## OpenMP
if FC_USE_OPENMP
AM_FCFLAGS += $(FCFLAGS_OPENMP)
AM_FFLAGS += $(FCFLAGS_OPENMP)
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
