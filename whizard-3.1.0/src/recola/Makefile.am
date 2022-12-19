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

## The files in this directory interface the Recola amplitude calculator

## We create a library which is still to be combined with auxiliary libs.
noinst_LTLIBRARIES = libwo_recola.la
check_LTLIBRARIES = libwo_recola_ut.la

if RECOLA_AVAILABLE
libwo_recola_la_SOURCES = \
  recola_wrapper.f90 \
  prc_recola.f90 \
  prc_recola_sub.f90
SOURCE_DUMMY_FILES = recola_wrapper_dummy.f90
else
libwo_recola_la_SOURCES = \
  recola_wrapper_dummy.f90 \
  prc_recola.f90 \
  prc_recola_sub.f90
SOURCE_DUMMY_FILES = recola_wrapper.f90
endif

libwo_recola_ut_la_SOURCES = \
  prc_recola_uti.f90 prc_recola_ut.f90

## Omitting this would exclude it from the distribution
dist_noinst_DATA = recola.nw

# Modules and installation
# Dump module names into file Modules
execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = \
  recola_wrapper.$(FCMOD) prc_recola.$(FCMOD)

libwo_recola_Modules = prc_recola recola_wrapper \
		    prc_recola_uti prc_recola_ut

if RECOLA_AVAILABLE
  prc_recola.lo prc_recola.$(FCMOD): recola_wrapper.lo
  recola_wrapper.$(FCMOD): recola_wrapper.lo
else
  prc_recola.lo prc_recola.$(FCMOD): recola_wrapper_dummy.lo
  recola_wrapper.$(FCMOD): recola_wrapper_dummy.lo
endif

Modules: Makefile
	@for module in $(libwo_recola_Modules); do \
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
  ../physics/Modules \
  ../qft/Modules \
  ../types/Modules \
  ../particles/Modules \
  ../matrix_elements/Modules \
  ../me_methods/Modules \
  ../variables/Modules \
  ../beams/Modules

$(module_lists):
	$(MAKE) -C `dirname $@` Modules

Module_dependencies.sed: $(module_lists)
	@rm -f $@
	echo 's/, *only:.*//' >> $@
	echo 's/, *&//' >> $@
	echo 's/, *.*=>.*//' >> $@
	echo 's/$$/.lo/' >> $@
	for list in $^; do \
	  dir="`dirname $$list`"; \
	  for mod in `cat $$list`; do \
	    echo 's!: '$$mod'.lo$$!': $$dir/$$mod'.lo!' >> $@; \
	  done ; \
	done

DISTCLEANFILES = Module_dependencies.sed


# The following line just says
#    include Makefile.depend
# but in a portable fashion (depending on automake's AM_MAKE_INCLUDE
@am__include@ @am__quote@Makefile.depend@am__quote@

Makefile.depend: Module_dependencies.sed
Makefile.depend: $(libwo_recola_la_SOURCES) ${libwo_recola_ut_la_SOURCES}
	@rm -f $@
	for src in $^; do \
	  module="`basename $$src | sed 's/\.f[90][0358]//'`"; \
	  grep '^ *use ' $$src \
	    | grep -v '!NODEP!' \
	    | sed -e 's/^ *use */'$$module'.lo: /' \
	          -f Module_dependencies.sed; \
	done > $@

DISTCLEANFILES += Makefile.depend 

### Disabled, explicit rule above
# Fortran90 module files are generated at the same time as object files
#.lo.$(FCMOD):#
#	@:
#	touch $@

AM_FCFLAGS = -I../basics -I../utilities -I../testing -I../system -I../combinatorics -I../parsing -I../physics -I../qft -I../expr_base -I../types -I../particles -I../matrix_elements -I../me_methods -I../variables -I../beams -I../lhapdf -I../pdf_builtin -I../fastjet 

########################################################################
prc_recola_sub.lo: prc_recola.lo

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

if RECOLA_AVAILABLE
AM_FCFLAGS += $(RECOLA_INCLUDES)
endif

########################################################################
## Non-standard targets and dependencies

## (Re)create F90 sources from NOWEB source.
if NOWEB_AVAILABLE

PRELUDE = $(top_srcdir)/src/noweb-frame/whizard-prelude.nw
POSTLUDE = $(top_srcdir)/src/noweb-frame/whizard-postlude.nw

recola.stamp: $(PRELUDE) $(srcdir)/recola.nw $(POSTLUDE)
	@rm -f recola.tmp
	@touch recola.tmp
	for src in $(libwo_recola_la_SOURCES) $(libwo_recola_ut_la_SOURCES) $(SOURCE_DUMMY_FILES); do \
	  $(NOTANGLE) -R[[$$src]] $^ | $(CPIF) $$src; \
        done
	@mv -f recola.tmp recola.stamp

$(libwo_recola_la_SOURCES) $(libwo_recola_ut_la_SOURCES) $(SOURCE_DUMMY_FILES): recola.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f recola.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) recola.stamp; \
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
	-rm -f recola.stamp recola.tmp
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
