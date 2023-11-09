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

bin_PROGRAMS = circe1_plot circe1_fit circe1_sample \
	circe1_minuit1 circe1_minuit2 circe1_int

bin_SCRIPTS = circe1_fit.run circe1_minuit1.run circe1_minuit2.run \
	circe1_int.run

KINDS=$(top_builddir)/circe1/src/kinds.lo

circe1_plot_SOURCES = circe1_plot.f90
circe1_plot_LDADD = $(KINDS) $(top_builddir)/circe1/src/libcirce1.la

circe1_sample_SOURCES = circe1_sample.f90
circe1_sample_LDADD = $(KINDS) $(top_builddir)/circe1/src/libcirce1.la

circe1_fit_SOURCES = circe1_fit.f90
circe1_fit_LDADD = $(KINDS) $(top_builddir)/circe1/minuit/libminuit.la \
	$(top_builddir)/circe1/src/libcirce1.la

circe1_int_SOURCES = circe1_int.f90
circe1_int_LDADD = $(KINDS) $(top_builddir)/circe1/src/libcirce1.la

circe1_minuit1_SOURCES = circe1_minuit1.f90
circe1_minuit1_LDADD = $(KINDS) $(top_builddir)/circe1/minuit/libminuit.la \
	$(top_builddir)/circe1/src/libcirce1.la

circe1_minuit2_SOURCES = circe1_minuit2.f90
circe1_minuit2_LDADD = $(KINDS) $(top_builddir)/circe1/minuit/libminuit.la \
	$(top_builddir)/circe1/src/libcirce1.la

EXTRA_DIST = $(circe1_plot_SOURCES) $(circe1_fit_SOURCES) \
	$(circe1_sample_SOURCES) $(circe1_int_SOURCES) \
	$(bin_SCRIPTS:.run=.sh)

SUFFIXES = .nw .lo .$(FCMOD) .sh. .run

.sh.run:
	@rm -f @$
	$(SED) 's|@name@|$*|g' $< >$@
	chmod +x $@

# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
	@:

AM_FFLAGS = 
AM_FCFLAGS = -I$(top_builddir)/circe1/src

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

BRACKET = $(top_srcdir)/circe1/src/prelude.nw $(top_srcdir)/circe1/src/postlude.nw
TRIPLE = $(top_srcdir)/circe1/src/prelude.nw $< $(top_srcdir)/circe1/src/postlude.nw

if NOWEB_AVAILABLE
circe1_plot.f90: $(top_srcdir)/circe1/src/circe1.nw $(BRACKET)
	cat $(TRIPLE) | $(NOTANGLE) -R'[[$@]]' > $@

circe1_sample.f90: $(top_srcdir)/circe1/src/circe1.nw $(BRACKET)
	cat $(TRIPLE) | $(NOTANGLE) -R'[[$@]]' > $@
endif

########################################################################

if NOWEB_AVAILABLE
circe1_fit.f90: $(top_srcdir)/circe1/src/minuit.nw $(BRACKET)
	cat $(TRIPLE) | $(NOTANGLE) -R'[[$@]]' > $@

circe1_fit.sh: $(top_srcdir)/circe1/src/minuit.nw $(BRACKET)
	cat $(TRIPLE) | $(NOTANGLE) -R'[[$@]]' > $@

circe1_minuit1.f90: $(top_srcdir)/circe1/src/minuit.nw $(BRACKET)
	cat $(TRIPLE) | $(NOTANGLE) -R'[[$@]]' > $@

circe1_minuit2.f90: $(top_srcdir)/circe1/src/minuit.nw $(BRACKET)
	cat $(TRIPLE) | $(NOTANGLE) -R'[[$@]]' > $@

circe1_minuit1.sh: $(top_srcdir)/circe1/src/minuit.nw $(BRACKET)
	cat $(TRIPLE) | $(NOTANGLE) -R'[[$@]]' > $@

circe1_minuit2.sh: $(top_srcdir)/circe1/src/minuit.nw $(BRACKET)
	cat $(TRIPLE) | $(NOTANGLE) -R'[[$@]]' > $@

circe1_int.f90: $(top_srcdir)/circe1/src/interpol.nw $(BRACKET)
	cat $(TRIPLE) | $(NOTANGLE) -R'[[$@]]' > $@

circe1_int.sh: $(top_srcdir)/circe1/src/interpol.nw $(BRACKET)
	cat $(TRIPLE) | $(NOTANGLE) -R'[[$@]]' > $@
endif

########################################################################
# The following line just says
#    include Makefile.depend
# but in a portable fashion (depending on automake's AM_MAKE_INCLUDE
########################################################################
@am__include@ @am__quote@Makefile.depend@am__quote@

Makefile.depend: $(circe1_plot_SOURCES) $(circe1_fit_SOURCES) \
	$(circe1_sample_SOURCES) $(circe1_int_SOURCES)
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

DISTCLEANFILES = Makefile.depend

########################################################################
## Non-standard cleanup tasks
## Remove sources that can be recreated using NOWEB
if NOWEB_AVAILABLE
maintainer-clean-noweb:
	-rm -f *.f90 *.c circe1_fit.sh circe1_minuit1.sh \
		circe1_minuit2.sh circe1_int.sh
endif
.PHONY: maintainer-clean-noweb

## Remove those sources also if builddir and srcdir are different
if NOWEB_AVAILABLE
clean-noweb:
	test "$(srcdir)" != "." && rm -f *.f90 *.c circe1_fit.sh \
	circe1_minuit1.sh circe1_minuit2.sh circe1_int.sh || true
endif
.PHONY: clean-noweb

## Remove F90 module files
clean-local: clean-noweb
	-rm -f *.run
	-rm -f *.$(FCMOD) *.g90 circemacs.mp4
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
	-test ! -f circe1_plot && $(MAKE) circe1_plot

###           module="`basename $$src | sed 's/\.f[90][0358]//'`"; 
########################################################################
# MPI
########################################################################
###
### # The -mismatch_all is for mpi_send() etc.
### MPIFC = mpif90
### MPIFCFLAGS = # -mismatch_all

