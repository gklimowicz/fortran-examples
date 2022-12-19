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

## The files in this directory interface the BLHA amplitude calculator

## We create a library which is still to be combined with auxiliary libs.
noinst_LTLIBRARIES = libblha.la
check_LTLIBRARIES = libblha_ut.la

COMMON_F90 = \
  blha_olp_interfaces.f90 \
  blha_config.f90
MPI_F90 = \
  blha_config_sub.f90_mpi
SERIAL_F90 = \
  blha_config_sub.f90_serial
BLHA_SUBMODULES = \
  blha_olp_interfaces_sub.f90

BLHA_MODULES = \
  $(COMMON_F90)

EXTRA_DIST = \
  $(COMMON_F90) \
  $(BLHA_SUBMODULES) \
  $(SERIAL_F90) \
  $(MPI_F90)

nodist_libblha_la_SOURCES = \
  $(BLHA_MODULES) \
  blha_config_sub.f90 \
  $(BLHA_SUBMODULES)

DISTCLEANFILES = blha_config_sub.f90

if FC_USE_MPI
blha_config_sub.f90: blha_config_sub.f90_mpi
	-cp -f $< $@
else
blha_config_sub.f90: blha_config_sub.f90_serial
	-cp -f $< $@
endif

libblha_ut_la_SOURCES = \
  blha_uti.f90 blha_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = blha.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  blha_olp_interfaces.$(FCMOD) \
  blha_config.$(FCMOD)

libblha_Modules = $(nodist_libblha_la_SOURCES:.f90=) $(libblha_ut_la_SOURCES:.f90=)
Modules: Makefile
	@for module in $(libblha_Modules); do \
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
  ../combinatorics/Modules \
  ../parsing/Modules \
  ../physics/Modules \
  ../qft/Modules \
  ../expr_base/Modules \
  ../types/Modules \
  ../variables/Modules \
  ../model_features/Modules \
  ../matrix_elements/Modules \
  ../particles/Modules \
  ../threshold/Modules \
  ../beams/Modules \
  ../me_methods/Modules

include_modules_bare = ${module_lists:/Modules=}
include_modules = ${include_modules_bare:../%=-I../%}

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(nodist_libblha_la_SOURCES) $(libblha_ut_la_SOURCES)
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

DISTCLEANFILES += Module_dependencies.sed

# The following line just says
#    include Makefile.depend
# but in a portable fashion (depending on automake's AM_MAKE_INCLUDE
@am__include@ @am__quote@Makefile.depend@am__quote@

Makefile.depend: Module_dependencies.sed
Makefile.depend: $(nodist_libblha_la_SOURCES) $(libblha_ut_la_SOURCES)
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

AM_FCFLAGS = $(include_modules) -I../fastjet -I../pdf_builtin -I../lhapdf

########################################################################
# For the moment, the submodule dependencies will be hard-coded
blha_config_sub.lo: blha_config.lo
blha_olp_interfaces_sub.lo: blha_olp_interfaces.lo

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

FILTER = -filter "sed 's/defn MPI:/defn/'"

PRELUDE = $(top_srcdir)/src/noweb-frame/whizard-prelude.nw
POSTLUDE = $(top_srcdir)/src/noweb-frame/whizard-postlude.nw

blha.stamp: $(PRELUDE) $(srcdir)/blha.nw $(POSTLUDE)
	@rm -f blha.tmp
	@touch blha.tmp
	for src in $(COMMON_F90) $(libblha_ut_la_SOURCES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	for src in $(BLHA_SUBMODULES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
	done
	for src in $(SERIAL_F90:.f90_serial=.f90); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src'_serial'; \
        done
	for src in $(MPI_F90:.f90_mpi=.f90); do \
	  $(NOTANGLE) -R[[$$src]] $(FILTER) $^ | $(CPIF) $$src'_mpi'; \
        done
	@mv -f blha.tmp blha.stamp

$(COMMON_F90) $(BLHA_SUBMODULES) $(SERIAL_F90) $(MPI_F90) $(libblha_ut_la_SOURCES): blha.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f blha.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) blha.stamp; \
	fi

endif


########################################################################
## Non-standard cleanup tasks
## Remove sources that can be recreated using NOWEB
if NOWEB_AVAILABLE
maintainer-clean-noweb:
	-rm -f *.f90 *.f90_serial *.f90_mpi *.c
endif
.PHONY: maintainer-clean-noweb

## Remove those sources also if builddir and srcdir are different
if NOWEB_AVAILABLE
clean-noweb:
	test "$(srcdir)" != "." && rm -f *.f90 *.f90_serial *.f90_mpi *.c || true
endif
.PHONY: clean-noweb

## Remove F90 module files
clean-local: clean-noweb
	-rm -f blha.stamp blha.tmp
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
