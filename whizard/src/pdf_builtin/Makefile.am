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

noinst_LTLIBRARIES = libpdf_builtin.la

PDF_MODULES = \
  cteq6pdf.f90 mrst2004qed.f90 mstwpdf.f90 ct10pdf.f90 \
  CJpdf.f90 ct14pdf.f90 ct18pdf.f90 pdf_builtin.f90

PDF_SUBMODULES = \
  pdf_builtin_sub.f90

libpdf_builtin_la_SOURCES = \
  $(PDF_MODULES) \
  $(PDF_SUBMODULES)

$(libpdf_builtin_la_OBJECTS): \
	../basics/kinds.$(FCMOD) \
	../basics/iso_varying_string.$(FCMOD) \
	../basics/io_units.$(FCMOD) \
	../utilities/format_utils.$(FCMOD) \
	../system/diagnostics.$(FCMOD)

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  pdf_builtin.$(FCMOD) \
  cteq6pdf.$(FCMOD) \
  cj_pdf.$(FCMOD) \
  mrst2004qed.$(FCMOD) \
  mstwpdf.$(FCMOD) \
  ct10pdf.$(FCMOD) \
  ct14pdf.$(FCMOD) \
  ct18pdf.$(FCMOD)

SUFFIXES = .lo .$(FCMOD)

# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
        @:
#       touch $@

AM_FCFLAGS = -I../basics -I../utilities -I../system

########################################################################
# For the moment, the submodule dependencies will be hard-coded
pdf_builtin_sub.lo: pdf_builtin.lo

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

pdf_builtin.lo: cteq6pdf.lo mrst2004qed.lo mstwpdf.lo ct10pdf.lo CJpdf.lo ct14pdf.lo ct18pdf.lo

pdf_builtin.$(FCMOD): pdf_builtin.lo
cteq6pdf.$(FCMOD): cteq6pdf.lo
cj_pdf.$(FCMOD): CJpdf.lo
mrst2004qed.$(FCMOD): mrst2004qed.lo
mstwpdf.$(FCMOD): mstwpdf.lo
ct10pdf.$(FCMOD): ct10pdf.lo
ct14pdf.$(FCMOD): ct14pdf.lo
ct18pdf.$(FCMOD): ct18pdf.lo

########################################################################
## Non-standard cleanup tasks

## Remove F90 module files
clean-local:
	-rm -f *.$(FCMOD)
if FC_SUBMODULES
	-rm -f *.smod *.sub
endif

## Remove backup files
maintainer-clean-local:
	-rm -f *~
