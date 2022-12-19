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

## The files in this directory implement the FKS subtraction scheme.

## We create a library which is still to be combined with auxiliary libs.
noinst_LTLIBRARIES = libfks.la
check_LTLIBRARIES = libfks_ut.la

libfks_la_SOURCES = \
  $(FKS_MODULES) \
  $(FKS_SUBMODULES)

FKS_MODULES = \
  fks_regions.f90 \
  nlo_data.f90 \
  virtual.f90 \
  real_subtraction.f90 \
  dglap_remnant.f90 \
  dispatch_fks.f90

FKS_SUBMODULES = \
  fks_regions_sub.f90 \
  nlo_data_sub.f90 \
  virtual_sub.f90 \
  real_subtraction_sub.f90 \
  dglap_remnant_sub.f90 \
  dispatch_fks_sub.f90

libfks_ut_la_SOURCES = \
  fks_regions_uti.f90 fks_regions_ut.f90 \
  real_subtraction_uti.f90 real_subtraction_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = fks.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ${FKS_MODULES:.f90=.$(FCMOD)}

libfks_Modules = ${FKS_MODULES:.f90=} ${libfks_ut_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(libfks_Modules); do \
          echo $$module >> $@.new; \
        done
	@if diff $@ $@.new -q >/dev/null; then \
          rm $@.new; \
	else \
          mv $@.new $@; echo "Modules updated"; \
	fi
BUILT_SOURCES = Modules

## Fortran module dependencies
# Get module lists from other directories
module_lists = \
  ../basics/Modules \
  ../utilities/Modules \
  ../testing/Modules \
  ../system/Modules \
  ../parsing/Modules \
  ../combinatorics/Modules \
  ../physics/Modules \
  ../qft/Modules \
  ../model_features/Modules \
  ../types/Modules \
  ../particles/Modules \
  ../threshold/Modules \
  ../matrix_elements/Modules \
  ../beams/Modules \
  ../me_methods/Modules \
  ../phase_space/Modules \
  ../variables/Modules \
  ../blha/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(libfks_la_SOURCES) $(libfks_ut_la_SOURCES)
Module_dependencies.sed: $(module_lists)
	@rm -f $@
	echo 's/, *only:.*//' >> $@
	echo 's/, *&//' >> $@
	echo 's/, *.*=>.*//' >> $@
	echo 's/$$/.lo/' >> $@
	for list in $(module_lists); do \
		dir="`dirname $$list`"; \
		for mod in `cat $$list`; do \
			echo 's!: '$$mod'.lo$$!': $$dir/$$mod'.lo!' >> $@; \
		done \
	done

DISTCLEANFILES = Module_dependencies.sed

# The following line just says
#    include Makefile.depend
# but in a portable fashion (depending on automake's AM_MAKE_INCLUDE
@am__include@ @am__quote@Makefile.depend@am__quote@

Makefile.depend: Module_dependencies.sed
Makefile.depend: $(libfks_la_SOURCES) ${libfks_ut_la_SOURCES}
	@rm -f $@
	for src in $^; do \
	  module="`basename $$src | sed 's/\.f[90][0358]//'`"; \
	  grep '^ *use ' $$src \
	    | grep -v '!NODEP!' \
	    | sed -e 's/^ *use */'$$module'.lo: /' \
	          -f Module_dependencies.sed; \
	done > $@

DISTCLEANFILES += Makefile.depend

# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
	@:
#	touch $@

AM_FCFLAGS = -I../basics -I../utilities -I../testing -I../system -I../combinatorics -I../parsing -I../rng -I../combinatorics -I../physics -I../qft -I../model_features -I../expr_base -I../types -I../particles -I../matrix_elements -I../beams -I../me_methods -I../phase_space -I../variables -I../blha -I../threshold -I../lhapdf -I../pdf_builtin -I../fastjet

########################################################################
fks_regions_sub.lo: fks_regions.lo
nlo_data_sub.lo: nlo_data.lo
virtual_sub.lo: virtual.lo
real_subtraction_sub.lo: real_subtraction.lo
dglap_remnant_sub.lo: dglap_remnant.lo
dispatch_fks_sub.lo: dispatch_fks.lo

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
## Non-standard targets and dependencies

## (Re)create F90 sources from NOWEB source.
if NOWEB_AVAILABLE

PRELUDE = $(top_srcdir)/src/noweb-frame/whizard-prelude.nw
POSTLUDE = $(top_srcdir)/src/noweb-frame/whizard-postlude.nw

fks.stamp: $(PRELUDE) $(srcdir)/fks.nw $(POSTLUDE)
	@rm -f fks.tmp
	@touch fks.tmp
	for src in $(libfks_la_SOURCES) $(libfks_ut_la_SOURCES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	@mv -f fks.tmp fks.stamp

$(libfks_la_SOURCES) $(libfks_ut_la_SOURCES): fks.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f fks.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) fks.stamp; \
	fi

endif


########################################################################
## Non-standard cleanup tasks
## Remove sources that can be recreated using NOWEB
if NOWEB_AVAILABLE
maintainer-clean-noweb:
	-rm -f *.f90 *.c
endif
.PHONY: maintainer-clean-noweb

## Remove those sources also if builddir and srcdir are different
if NOWEB_AVAILABLE
clean-noweb:
	test "$(srcdir)" != "." && rm -f *.f90 *.c || true
endif
.PHONY: clean-noweb

## Remove F90 module files
clean-local: clean-noweb
	-rm -f fks.stamp fks.tmp
	-rm -f *.$(FCMOD)
if FC_SUBMODULES
	-rm -f *.smod *.sub
endif

## Remove backup files
maintainer-clean-backup:
	-rm -f *~
.PHONY: maintainer-clean-backup

## Register additional clean targets
maintainer-clean-local: maintainer-clean-noweb maintainer-clean-backup
