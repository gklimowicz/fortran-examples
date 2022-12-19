## Makefile.am -- Makefile for model-parameter modules in WHIZARD
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

## Build the threeshl bundle required for the 3SHL models.

noinst_LTLIBRARIES = libthreeshl_bundle.la

## We need both in the distribution

EXTRA_DIST = threeshl_bundle.f90 threeshl_bundle_lt.f90

libthreeshl_bundle_la_SOURCES = threeshl_bundle.f90

## Unfortunately, the NAG and intel compilers require these to be happy.
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
	tdefs.$(FCMOD) \
	threeshl.$(FCMOD) \
	tscript.$(FCMOD) \
	tglue.$(FCMOD)

# Fortran90 module files are generated at the same time as object files
$(nodist_execmod_HEADERS): threeshl_bundle.lo
	@:

########################################################################
## Dependencies on global WHIZARD modules

basicsdir = $(top_builddir)/src/basics

KINDS_MOD = kinds.$(FCMOD)
CONSTANTS_MOD = constants.$(FCMOD)

AM_FCFLAGS = -I$(basicsdir) 

$(libthreeshl_bundle_la_OBJECTS): \
  $(basicsdir)/$(KINDS_MOD) \
  $(basicsdir)/$(CONSTANTS_MOD)

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
clean-local:
	-rm -f *.$(FCMOD)
if FC_SUBMODULES
	-rm -f *.smod
endif

## Remove backup files
maintainer-clean-local:
	-rm -f *~
