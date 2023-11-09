# Makefile.am --
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

SOURCE_FILES = \
  circe1.f90

MODULE_FILES = \
  kinds.$(FCMOD) \
  circe1.$(FCMOD)

NOWEB_FILES = \
    prelude.nw postlude.nw circe1.nw \
    minuit.nw interpol.nw taorng.nw

EXTRA_DIST = $(NOWEB_FILES) $(SOURCE_FILES)

# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
	@:

lib_LTLIBRARIES = libcirce1.la
libcirce1_la_SOURCES = $(SOURCE_FILES)

AM_FFLAGS = 
AM_FCFLAGS = 

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
# noweb
########################################################################

BRACKET = prelude.nw postlude.nw
TRIPLE = $(srcdir)/prelude.nw $< $(srcdir)/postlude.nw
NOTANGLE_IT = \
  cat $(TRIPLE) | $(NOTANGLE) -R'[[$@]]' > $@

SUFFIXES = .nw .lo .$(FCMOD)

if NOWEB_AVAILABLE
.nw.f90:
	$(NOTANGLE_IT)
endif NOWEB_AVAILABLE

execmodcirce1dir = $(fmoddir)/circe1
nodist_execmodcirce1_HEADERS = $(MODULE_FILES)

execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = circe1.$(FCMOD)

########################################################################
# The following line just says
#    include Makefile.depend
# but in a portable fashion (depending on automake's AM_MAKE_INCLUDE
########################################################################
@am__include@ @am__quote@Makefile.depend@am__quote@

Makefile.depend: $(SOURCE_FILES) 
	@rm -f $@
	for src in $^; do \
          module="`basename $$src | sed 's/\.f90//'`"; \
	  grep '^ *use ' $$src \
	    | grep -v '!NODEP!' \
	    | sed -e 's/^ *use */'$$module'.lo: /' \
	          -e 's/, *only:.*//' \
	          -e 's/, *&//' \
	          -e 's/, *.*=>.*//' \
	          -e 's/ *$$/.lo/' \
				 $(depend_filter_extra) ; \
	done > $@

DISTCLEANFILES = Makefile.depend kinds.f90

########################################################################
## Non-standard cleanup tasks
## Remove sources that can be recreated using NOWEB

if NOWEB_AVAILABLE
maintainer-clean-noweb:
	-rm -f $(SOURCE_FILES)
endif
.PHONY: maintainer-clean-noweb

## Remove those sources also if builddir and srcdir are different
if NOWEB_AVAILABLE
clean-noweb:
	test "$(srcdir)" != "." && rm -f $(SOURCE_FILES) || true
endif
.PHONY: clean-noweb

## Remove F90 module files
clean-local: clean-noweb
	-rm -f *.$(FCMOD) *.g90
if FC_SUBMODULES
	-rm -f *.smod
endif

## Remove backup files
maintainer-clean-backup:
	-rm -f *~
.PHONY: maintainer-clean-backup

## Register additional clean targets
maintainer-clean-local: maintainer-clean-noweb maintainer-clean-backup

dist-hook:
	-test ! -f libcirce1.la && $(MAKE) libcirce1.la

###           module="`basename $$src | sed 's/\.f[90][0358]//'`"; 
########################################################################
# MPI
########################################################################
###
### # The -mismatch_all is for mpi_send() etc.
### MPIFC = mpif90
### MPIFCFLAGS = # -mismatch_all

