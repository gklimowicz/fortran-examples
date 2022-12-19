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

## The files in this directory implement physics definitions and functions
## for use in the WHIZARD generator.

## We create a library which is still to be combined with auxiliary libs.
noinst_LTLIBRARIES = libphysics.la
check_LTLIBRARIES = libphysics_ut.la

libphysics_la_SOURCES = \
  $(PHYSICS_MODULES) \
  $(PHYSICS_SUBMODULES)

PHYSICS_MODULES = \
  physics_defs.f90 \
  c_particles.f90 \
  lorentz.f90 \
  phs_points.f90 \
  sm_physics.f90 \
  sm_qcd.f90 \
  sm_qed.f90 \
  shower_algorithms.f90

PHYSICS_SUBMODULES = \
  physics_defs_sub.f90 \
  c_particles_sub.f90 \
  lorentz_sub.f90 \
  phs_points_sub.f90 \
  sm_physics_sub.f90 \
  sm_qcd_sub.f90 \
  sm_qed_sub.f90 \
  shower_algorithms_sub.f90

libphysics_ut_la_SOURCES = \
  sm_physics_uti.f90 sm_physics_ut.f90 \
  sm_qcd_uti.f90 sm_qcd_ut.f90 \
  sm_qed_uti.f90 sm_qed_ut.f90 \
  lorentz_uti.f90 lorentz_ut.f90 \
  phs_points_uti.f90 phs_points_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = physics.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  ${PHYSICS_MODULES:.f90=.$(FCMOD)}

# Submodules must not be included here
libphysics_Modules = \
  ${PHYSICS_MODULES:.f90=} \
  ${libphysics_ut_la_SOURCES:.f90=}
Modules: Makefile
	@for module in $(libphysics_Modules); do \
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
  ../combinatorics/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(libphysics_la_SOURCES) $(libphysics_ut_la_SOURCES)
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
Makefile.depend: $(libphysics_la_SOURCES) $(libphysics_ut_la_SOURCES)
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

AM_FCFLAGS = -I../basics -I../utilities -I../testing -I../system -I../combinatorics

########################################################################
# For the moment, the submodule dependencies will be hard-coded
physics_defs_sub.lo: physics_defs.lo
c_particles_sub.lo: c_particles.lo
lorentz_sub.lo: lorentz.lo
phs_points_sub.lo: phs_points.lo
sm_physics_sub.lo: sm_physics.lo
sm_qcd_sub.lo: sm_qcd.lo
sm_qed_sub.lo: sm_qed.lo
shower_algorithms_sub.lo: shower_algorithms.lo

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

physics.stamp: $(PRELUDE) $(srcdir)/physics.nw $(POSTLUDE)
	@rm -f physics.tmp
	@touch physics.tmp
	for src in $(libphysics_la_SOURCES) $(libphysics_ut_la_SOURCES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	@mv -f physics.tmp physics.stamp

$(libphysics_la_SOURCES) $(libphysics_ut_la_SOURCES): physics.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f physics.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) physics.stamp; \
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
	-rm -f physics.stamp physics.tmp
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
