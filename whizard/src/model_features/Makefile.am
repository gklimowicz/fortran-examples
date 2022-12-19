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

## The files in this directory handle process definitions and models

## We create a library which is still to be combined with auxiliary libs.
noinst_LTLIBRARIES = libmodel_features.la
check_LTLIBRARIES = libmodel_features_ut.la

libmodel_features_la_SOURCES = \
  $(MODEL_FEATURES_MODULES) \
  $(MODEL_FEATURES_SUBMODULES)

MODEL_FEATURES_MODULES = \
  auto_components.f90 \
  radiation_generator.f90 \
  eval_trees.f90 \
  models.f90 \
  slha_interface.f90

MODEL_FEATURES_SUBMODULES = \
  auto_components_sub.f90 \
  radiation_generator_sub.f90 \
  eval_trees_sub.f90 \
  models_sub.f90 \
  slha_interface_sub.f90

libmodel_features_ut_la_SOURCES = \
  auto_components_uti.f90 auto_components_ut.f90 \
  radiation_generator_uti.f90 radiation_generator_ut.f90 \
  eval_trees_uti.f90 eval_trees_ut.f90 \
  models_uti.f90 models_ut.f90 \
  slha_interface_uti.f90 slha_interface_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = model_features.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ${MODEL_FEATURES_MODULES:.f90=.$(FCMOD)}

# This must not include the submodules
libmodel_features_Modules = \
  ${MODEL_FEATURES_MODULES:.f90=} \
  ${libmodel_features_ut_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(libmodel_features_Modules); do \
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
  ../particles/Modules \
  ../variables/Modules \
  ../threshold/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(libmodel_features_la_SOURCES) \
  $(libmodel_features_ut_la_SOURCES)
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
Makefile.depend: $(libmodel_features_la_SOURCES) \
  $(libmodel_features_ut_la_SOURCES)
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

AM_FCFLAGS = -I../basics -I../utilities -I../testing -I../system -I../combinatorics -I../parsing -I../physics -I../qft -I../expr_base -I../types -I../particles -I../threshold -I../fastjet -I../variables

########################################################################
# For the moment, the submodule dependencies will be hard-coded
auto_components_sub.lo: auto_components.lo
radiation_generator_sub.lo: radiation_generator.lo
eval_trees_sub.lo: eval_trees.lo
models_sub.lo: models.lo
slha_interface_sub.lo: slha_interface.lo

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

model_features.stamp: $(PRELUDE) $(srcdir)/model_features.nw $(POSTLUDE)
	@rm -f model_features.tmp
	@touch model_features.tmp
	for src in \
          $(libmodel_features_la_SOURCES) \
          $(libmodel_features_ut_la_SOURCES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	@mv -f model_features.tmp model_features.stamp

$(libmodel_features_la_SOURCES) $(libmodel_features_ut_la_SOURCES): model_features.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f model_features.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) model_features.stamp; \
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
	-rm -f model_features.stamp model_features.tmp
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
