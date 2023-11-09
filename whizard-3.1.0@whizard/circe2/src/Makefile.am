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

# backwards compatibility
OCAML_BIGARRAY_COMPAT=@OCAML_BIGARRAY_COMPAT@
OCAML_BIGARRAY_CMA=@OCAML_BIGARRAY_CMA@
OCAML_BIGARRAY_CMXA=@OCAML_BIGARRAY_CMXA@

########################################################################

lib_LTLIBRARIES =
SOURCE_FILES =
MODULE_FILES =
EXTRA_SOURCE_FILES = bigarray_library.ml bigarray_module.ml bigarray_compat.mli
NOWEB_FILES = prelude.nw postlude.nw

########################################################################
# The CIRCE2 library proper
########################################################################

lib_LTLIBRARIES += libcirce2.la
libcirce2_la_SOURCES = circe2.f90
MODULE_FILES += circe2.$(FCMOD)
SOURCE_FILES += $(libcirce2_la_SOURCES)
NOWEB_FILES += circe2.nw

########################################################################
# Required for standalone compilation,
# otherwise provided by VAMP and/or WHIZARD
########################################################################

EXTRA_SOURCE_FILES +=
MODULE_FILES += kinds.$(FCMOD)

lib_LTLIBRARIES += libtaorng.la
libtaorng_la_SOURCES = tao_random_numbers.f90
MODULE_FILES += tao_random_numbers.$(FCMOD)
SOURCE_FILES += $(libtaorng_la_SOURCES)

lib_LTLIBRARIES += libtaorng_objs.la
libtaorng_objs_la_SOURCES = tao_random_objects.f90
MODULE_FILES += tao_random_objects.$(FCMOD)
SOURCE_FILES += $(libtaorng_objs_la_SOURCES)

########################################################################
# Install generated .mod files
#   once in 'circe2' and once in 'whizard' (main only)
########################################################################
execmodcircedir = $(fmoddir)/circe2
nodist_execmodcirce_HEADERS = $(MODULE_FILES)

execmoddir = $(fmoddir)/whizard
nodist_execmod_HEADERS = circe2.$(FCMOD)

########################################################################
# Testing and tools
########################################################################

bin_PROGRAMS = circe2_moments circe2_ls circe2_generate

circe2_moments_SOURCES = circe2_moments.f90
circe2_moments_LDADD = kinds.lo libcirce2.la libtaorng_objs.la libtaorng.la 
circe2_moments.o: $(MODULE_FILES)

circe2_ls_SOURCES = circe2_ls.f90
circe2_ls_LDADD = kinds.lo libcirce2.la libtaorng.la
circe2_ls.o: $(MODULE_FILES)

circe2_generate_SOURCES = circe2_generate.f90
circe2_generate_LDADD = kinds.lo libcirce2.la libtaorng_objs.la libtaorng.la 
circe2_generate.o: $(MODULE_FILES)

bin_SCRIPTS =

if OCAML_AVAILABLE
bin_SCRIPTS += $(CIRCE2_NATIVE)
endif OCAML_AVAILABLE

$(CIRCE2_NATIVE): $(CIRCE2_CMX)
$(CIRCE2_BYTECODE): $(CIRCE2_CMO)

if OCAML_AVAILABLE
all-local: $(CIRCE2_CMX) $(CIRCE2TOOL_CMX)
bytecode: $(CIRCE2_CMO) $(CIRCE2TOOL_CMO)
else
all-local:
bytecode:
endif

include $(top_srcdir)/circe2/src/Makefile.ocaml
include $(top_srcdir)/circe2/src/Makefile.sources

EXTRA_DIST = $(NOWEB_FILES) $(SOURCE_FILES) $(CIRCE2_CAML) \
	$(EXTRA_SOURCE_FILES) $(circe2_moments_SOURCES) \
	$(circe2_ls_SOURCES)  $(circe2_generate_SOURCES)

MYPRECIOUS = $(CIRCE2_DERIVED)

# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
	@:
#	touch $@

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

TRIPLE = $(srcdir)/prelude.nw $< $(srcdir)/postlude.nw
WEBS = $(srcdir)/prelude.nw $(srcdir)/circe2.nw $(srcdir)/postlude.nw
NOTANGLE_IT = \
  cat $(TRIPLE) | $(NOTANGLE) -R'[[$@]]' > $@

SUFFIXES += .nw .$(FCMOD)

if NOWEB_AVAILABLE
.nw.f90:
	$(NOTANGLE_IT)

circe2_ls.f90: circe2.nw
	cat $(WEBS) | $(NOTANGLE) -R'[[$@]]' > $@

circe2_generate.f90: circe2.nw
	cat $(WEBS) | $(NOTANGLE) -R'[[$@]]' > $@

circe2_moments.f90: circe2.nw
	cat $(WEBS) | $(NOTANGLE) -R'[[$@]]' > $@

tao_random_objects.f90: circe2.nw
	cat $(WEBS) | $(NOTANGLE) -R'[[$@]]' > $@

endif NOWEB_AVAILABLE

########################################################################
# O'Caml
########################################################################

if OCAML_AVAILABLE

bigarray_compat.ml: $(OCAML_BIGARRAY_COMPAT).ml
	cp -f $< $@

events.cmx: bigarray_compat.cmi bigarray_compat.cmx
events.cmo: bigarray_compat.cmi bigarray_compat.cmo

circe2.top: $(CIRCE2_CMO)
	$(OCAMLMKTOP) $(OCAMLFLAGS) -o $@ \
		unix.cma $(OCAML_BIGARRAY_CMA) $(CIRCE2_CMO)

lexer.mli: lexer.ml parser.cmi
	$(OCAMLC) -i $< | $(GREP) 'val token' >$@

events_lexer.mli: events_lexer.ml
	$(OCAMLC) -i $< | $(GREP) 'val token' >$@

endif OCAML_AVAILABLE


########################################################################
# The following line just says
#    include Makefile.depend
# but in a portable fashion (depending on automake's AM_MAKE_INCLUDE
########################################################################
@am__include@ @am__quote@Makefile.depend@am__quote@

Makefile.depend: $(SOURCE_FILES) $(circe2_moments_SOURCES) $(circe2_ls_SOURCES)  $(circe2_generate_SOURCES)
	@rm -f $@
	for src in $^; do \
          module="`basename $$src | sed 's/\.f90//'`"; \
	  grep '^ *use ' $$src \
	    | grep -v '!NODEP!' \
	    | sed -e 's/^ *use */'$$module'.lo: /' \
	          -e 's/, *only:.*//' \
	          -e 's/, *&//' \
	          -e 's/, *.*=>.*//' \
	          -e 's/ *$$/.lo/'; \
	  grep '^ *use ' $$src \
	    | grep -v '!NODEP!' \
	    | sed -e 's/^ *use */'$$module'.lo: /' \
	          -e 's/, *only:.*//' \
	          -e 's/, *&//' \
	          -e 's/, *.*=>.*//' \
	          -e 's/ *$$/.$$(FCMOD)/'; \
	done > $@

DISTCLEANFILES = Makefile.depend kinds.f90

if OCAML_AVAILABLE

@am__include@ @am__quote@Makefile.depend_ocaml@am__quote@

# echo lexer.mli: lexer.ml >>$@
Makefile.depend_ocaml: $(CIRCE2_SRC) $(CIRCE2TOOL_SRC)
	@if $(AM_V_P); then :; else echo "  OCAMLDEP " $@; fi
	@rm -f $@
	$(AM_V_at)$(OCAMLDEP) -I $(srcdir) $^ | sed 's,[^ 	]*/,,g' >$@
	echo parser.mli: parser.ml >>$@
	echo lexer.cmi: parser.cmi >>$@
	echo parser.cmi: syntax.cmi >>$@
	echo commands.cmi: parser.cmi lexer.cmi >>$@
	echo commands.cmo: parser.cmi lexer.cmi >>$@
	echo commands.cmx: parser.cmx lexer.cmx >>$@
	echo lexer.cmo: lexer.cmi >>$@
	echo lexer.cmx: lexer.cmi parser.cmx >>$@
	echo parser.cmo: parser.cmi syntax.cmi >>$@
	echo parser.cmx: parser.cmi syntax.cmi syntax.cmx >>$@
	echo events_lexer.cmo: events_lexer.cmi >>$@
	echo events_lexer.cmx: events_lexer.cmi >>$@

DISTCLEANFILES += Makefile.depend_ocaml

DISTCLEANFILES += $(CIRCE2_DERIVED)

DISTCLEANFILES += bigarray_compat.ml

endif OCAML_AVAILABLE

########################################################################
## Non-standard cleanup tasks
## Remove sources that can be recreated using NOWEB
.PRECIOUS = $(MYPRECIOUS)

if NOWEB_AVAILABLE
maintainer-clean-noweb:
	-rm -f $(SOURCE_FILES) $(circe2_moments_SOURCES) $(circe2_ls_SOURCES)  $(circe2_generate_SOURCES)
endif
.PHONY: maintainer-clean-noweb

## Remove those sources also if builddir and srcdir are different
if NOWEB_AVAILABLE
clean-noweb:
	test "$(srcdir)" != "." && rm -f $(SOURCE_FILES) $(circe2_moments_SOURCES) $(circe2_ls_SOURCES)  $(circe2_generate_SOURCES) || true
endif
.PHONY: clean-noweb

## Remove F90 module files
clean-local: clean-noweb
	-rm -f *.cm[aiox] *.cmxa *.[ao] *.l[oa] *.$(FCMOD) \
	*.g90 $(CIRCE2_NATIVE) $(CIRCE2_BYTECODE) $(CIRCE2_DERIVED)
if FC_SUBMODULES
	-rm -f *.smod
endif

## Remove backup files
maintainer-clean-backup:
	-rm -f *~
.PHONY: maintainer-clean-backup

## Register additional clean targets
maintainer-clean-local: maintainer-clean-noweb maintainer-clean-backup

###           module="`basename $$src | sed 's/\.f[90][0358]//'`"; 
########################################################################
# MPI
########################################################################
###
### # The -mismatch_all is for mpi_send() etc.
### MPIFC = mpif90
### MPIFCFLAGS = # -mismatch_all

