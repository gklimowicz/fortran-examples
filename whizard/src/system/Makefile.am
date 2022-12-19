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

## The files in this directory implement OS interactions of WHIZARD

## We create a library which is still to be combined with auxiliary libs.
noinst_LTLIBRARIES = libsystem.la
check_LTLIBRARIES = libsystem_ut.la

COMMON_F90 = \
  system_defs.f90 \
  signal_interface.c \
  sprintf_interface.c \
  diagnostics.f90 \
  os_interface.f90 \
  formats.f90 \
  cputime.f90
MPI_F90 = \
  os_interface_sub.f90_mpi
SERIAL_F90 = \
  os_interface_sub.f90_serial
SYSTEM_SUBMODULES = \
  diagnostics_sub.f90 \
  formats_sub.f90 \
  cputime_sub.f90

SYSTEM_MODULES = \
  system_dependencies.f90 \
  debug_master.f90 \
  $(COMMON_F90)

EXTRA_DIST = \
  $(COMMON_F90) \
  $(SYSTEM_SUBMODULES) \
  $(SERIAL_F90) \
  $(MPI_F90)

nodist_libsystem_la_SOURCES = \
  $(SYSTEM_MODULES) \
  $(SYSTEM_SUBMODULES) \
  os_interface_sub.f90

DISTCLEANFILES = os_interface_sub.f90

if FC_USE_MPI
os_interface_sub.f90: os_interface_sub.f90_mpi
	-cp -f $< $@
else
os_interface_sub.f90: os_interface_sub.f90_serial
	-cp -f $< $@
endif

libsystem_ut_la_SOURCES = \
  os_interface_uti.f90 os_interface_ut.f90 \
  formats_uti.f90 formats_ut.f90 \
  cputime_uti.f90 cputime_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = system.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  system_dependencies.$(FCMOD) \
  debug_master.$(FCMOD) \
  system_defs.$(FCMOD) \
  cputime.$(FCMOD) \
  diagnostics.$(FCMOD) \
  formats.$(FCMOD) \
  os_interface.$(FCMOD)

libsystem_Modules = \
  $(SYSTEM_MODULES:.f90=) \
  $(libsystem_ut_la_SOURCES:.f90=)
Modules: Makefile
	@for module in $(libsystem_Modules); do \
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
  ../testing/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(nodist_libsystem_la_SOURCES) $(libsystem_ut_la_SOURCES)
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
Makefile.depend: $(nodist_libsystem_la_SOURCES) $(libsystem_ut_la_SOURCES)
	@rm -f $@
	for src in $^; do \
	  module="`basename $$src | sed 's/\.f[90][0358]//'`"; \
	  grep '^ *use ' $$src \
	    | grep -v '!NODEP!' \
	    | sed -e 's/^ *use */'$$module'.lo: /' \
	          -f Module_dependencies.sed; \
	done > $@

DISTCLEANFILES += Makefile.depend

SUFFIXES = .lo .$(FCMOD)

# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
	@:
#	touch $@

AM_FCFLAGS = -I../basics -I../utilities -I../testing

########################################################################
# For the moment, the submodule dependencies will be hard-coded
diagnostics_sub.lo: diagnostics.lo
os_interface_sub.lo: os_interface.lo
formats_sub.lo: formats.lo
cputime_sub.lo: cputime.lo

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

## MPI
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

system.stamp: $(PRELUDE) $(srcdir)/system.nw $(POSTLUDE)
	@rm -f system.tmp
	@touch system.tmp
	for src in $(COMMON_F90) $(libsystem_ut_la_SOURCES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	for src in $(SYSTEM_SUBMODULES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	for src in $(SERIAL_F90:.f90_serial=.f90); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src'_serial'; \
        done
	for src in $(MPI_F90:.f90_mpi=.f90); do \
	  $(NOTANGLE) -R[[$$src]] $(FILTER) $^ | $(CPIF) $$src'_mpi'; \
        done
	@mv -f system.tmp system.stamp

$(COMMON_F90) $(SYSTEM_SUBMODULES) $(SERIAL_F90) $(MPI_F90) $(libsystem_ut_la_SOURCES): system.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f system.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) system.stamp; \
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
	-rm -f system.stamp system.tmp
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
