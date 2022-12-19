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

noinst_LTLIBRARIES = libwo_stdhep.la

libwo_stdhep_la_SOURCES = \
   mcf_hepev4_xdr.c \
   mcf_hepup_xdr.c \
   mcf_stdcm1_xdr.c \
   mcf_Stdhep_xdr.c \
   stdhep_mcfio.c \
   stdhep_internal_utils.c \
   stdtcopy.f \
   stdversn.f \
   stdxend.f \
   stdxrd.f \
   stdxrinit.f \
   stdxropen.f \
   stdxwcm1.f \
   stdxwevt.f \
   stdxwevtlh.f \
   stdxwevtup.f \
   stdxwinit.f \
   stdxwopen.f \
   stdxwrt.f \
   stdzero.f \
   hepev4.h \
   hepeup.h \
   heprup.h \
   stdcm1.h \
   stdcnt.h \
   stdhd.h \
   stdhep.h \
   stdhep_mcfio.h \
   stdlun.h \
   stdtmp.h \
   stdver.h \
   hepev4.inc \
   hepeup.inc \
   heprup.inc \
   mcfio.inc \
   stdcm1.inc \
   stdcnt.inc \
   stdhd.inc \
   stdhep.inc \
   stdlun.inc \
   stdtmp.inc \
   stdver.inc

## The include files are not automatically found (by certain Fortran compilers)
## RPC flags need to be included if SunRPC is absent
AM_FCFLAGS = -I$(top_srcdir)/stdhep
AM_FFLAGS = -I$(top_srcdir)/stdhep
AM_CFLAGS = $(RPC_CFLAGS) -I$(top_srcdir)/mcfio 

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
