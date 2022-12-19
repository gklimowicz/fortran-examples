# Makefile.am --
########################################################################
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

depend_filter_extra = -e 's/iso_fortran_env\.lo//'

NOWEB_SOURCE_FILES = \
  divisions.f90 \
  vamp.f90 \
  exceptions.f90 \
  tao_random_numbers.f90 \
  tao52_random_numbers.f90 \
  specfun.f90 \
  vamp_stat.f90 \
  histograms.f90 \
  utils.f90 \
  linalg.f90 \
  products.f90 \
  $(sources_extra)

SOURCE_FILES = $(NOWEB_SOURCE_FILES)

MODULE_FILES = \
  divisions.$(FCMOD) \
  vamp.$(FCMOD) \
  kinds.$(FCMOD) \
  exceptions.$(FCMOD) \
  tao_random_numbers.$(FCMOD) \
  tao52_random_numbers.$(FCMOD) \
  specfun.$(FCMOD) \
  vamp_stat.$(FCMOD) \
  histograms.$(FCMOD) \
  utils.$(FCMOD) \
  linalg.$(FCMOD) \
  products.$(FCMOD) \
  $(modules_extra) \
  vamp_grid_type.$(FCMOD) \
  vamp_rest.$(FCMOD) \
  vamp_equivalences.$(FCMOD)

NOWEB_FILES = \
    prelude.nw postlude.nw \
    divisions.nw vamp.nw \
    vamp_kinds.nw constants.nw exceptions.nw \
    tao_random_numbers.nw \
    specfun.nw vamp_stat.nw histograms.nw utils.nw \
    linalg.nw products.nw kinematics.nw coordinates.nw \
    vampi.nw mpi90.nw vamp_test.nw vamp_test0.nw \
    application.nw sample.nw

# Modules and installation
# Install twice: once for vamp
execvampmoddir = $(fmoddir)/vamp
nodist_execvampmod_HEADERS = $(MODULE_FILES)

# Once for whizard
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = $(MODULE_FILES)

lib_LTLIBRARIES = libvamp.la
libvamp_la_SOURCES = $(SOURCE_FILES)

EXTRA_DIST = $(NOWEB_FILES) $(SOURCE_FILES) iso_fortran_env_stub.f90

# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
	@:
#	touch $@

# Explicit dependencies
vamp_grid_type.$(FCMOD): vamp.lo
vamp_rest.$(FCMOD): vamp.lo
vamp_equivalences.$(FCMOD): vamp.lo

########################################################################
# Fortran Compilers
########################################################################

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

########################################################################
# noweb
########################################################################

NOTANGLE_IT = \
  $(NOTANGLE) -R'[[$(@F)]]' $(srcdir)/prelude.nw $^ | $(CPIF) $@

SUFFIXES = .nw .lo .$(FCMOD) 

.nw.f90:
	$(NOTANGLE_IT)

tao52_random_numbers.f90: tao_random_numbers.nw
	$(NOTANGLE_IT)

########################################################################
# The following line just says
#    include Makefile.depend
# but in a portable fashion (depending on automake's AM_MAKE_INCLUDE
########################################################################
@am__include@ @am__quote@Makefile.depend@am__quote@

Makefile.depend: $(SOURCE_FILES)
	@rm -f $@
	for src in $^; do \
	  module="`basename $$src | sed 's/\.f[90][0358]//'`"; \
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
	-rm -f $(NOWEB_SOURCE_FILES)
endif
.PHONY: maintainer-clean-noweb

## Remove those sources also if builddir and srcdir are different
if NOWEB_AVAILABLE
clean-noweb:
	test "$(srcdir)" != "." && rm -f $(NOWEB_SOURCE_FILES) || true
endif
.PHONY: clean-noweb

## Remove F90 module files
clean-local: clean-noweb
	-rm -f *.$(FCMOD)
if FC_SUBMODULES
	-rm -f *.smod
endif

## Remove backup files
maintainer-clean-backup:
	-rm -f *~
.PHONY: maintainer-clean-backup

## Register additional clean targets
maintainer-clean-local: maintainer-clean-noweb maintainer-clean-backup

########################################################################
# MPI
########################################################################
###
### # The -mismatch_all is for mpi_send() etc.
### MPIFC = mpif90
### MPIFCFLAGS = # -mismatch_all
