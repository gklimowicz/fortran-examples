# Makefile.am -- Makefile for O'Mega within and without WHIZARD
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

# Build the O'Mega Fortran library using libtool
# (?use pkglib_ instead of lib_ to make the -rpath and *.lai business work ...)
lib_LTLIBRARIES = libomega_core.la
execmoddir = $(fmoddir)/omega
nodist_execmod_HEADERS = $(OMEGALIB_MOD)

libomega_core_la_SOURCES = $(OMEGALIB_F90)

EXTRA_DIST = \
    $(OMEGA_CONFIG_MLI) $(OMEGA_CAML) \
    omegalib.nw $(OMEGALIB_F90)

OMEGA_CMXA = omega_core.cmxa omega_targets.cmxa omega_models.cmxa
OMEGA_CMA = $(OMEGA_CMXA:.cmxa=.cma)

if OCAML_AVAILABLE 
all-local: $(OMEGA_CMXA) $(OMEGA_APPLICATIONS_CMX)
bytecode: $(OMEGA_CMA) $(OMEGA_APPLICATIONS_CMO)
else
all-local: 
bytecode: 
endif

# Compiled interfaces and libraries for out-of-tree compilation of models
if OCAML_AVAILABLE
camllibdir = $(libdir)/omega/caml
nodist_camllib_DATA = \
	omega.cmi fusion.cmi targets.cmi coupling.cmi modeltools.cmi color.cmi \
	options.cmi model.cmi \
	omega_core.cmxa omega_core.a omega_targets.cmxa omega_targets.a \
	charges.cmi
endif

########################################################################

include $(top_srcdir)/omega/src/Makefile.ocaml
include $(top_srcdir)/omega/src/Makefile.sources

if OCAML_AVAILABLE
config.cmo config.cmx: config.cmi

omega_core.a: omega_core.cmxa
omega_core.cmxa: $(OMEGA_CORE_CMX)
	@if $(AM_V_P); then :; else echo "  OCAMLOPT " $@; fi
	$(AM_V_at)$(OCAMLOPT) $(OCAMLFLAGS) $(OCAMLOPTFLAGS) -a -o $@ $^

omega_core.cma: $(OMEGA_CORE_CMO)
	@if $(AM_V_P); then :; else echo "  OCAMLC   " $@; fi
	$(AM_V_at)$(OCAMLC) $(OCAMLFLAGS) -a -o $@ $^

omega_targets.a: omega_targets.cmxa
omega_targets.cmxa: $(OMEGA_TARGETS_CMX)
	@if $(AM_V_P); then :; else echo "  OCAMLOPT " $@; fi
	$(AM_V_at)$(OCAMLOPT) $(OCAMLFLAGS) $(OCAMLOPTFLAGS) -a -o $@ $^

omega_targets.cma: $(OMEGA_TARGETS_CMO)
	@if $(AM_V_P); then :; else echo "  OCAMLC   " $@; fi
	$(AM_V_at)$(OCAMLC) $(OCAMLFLAGS) -a -o $@ $^

omega_models.cmxa: $(OMEGA_MODELS_CMX)
	@if $(AM_V_P); then :; else echo "  OCAMLOPT " $@; fi
	$(AM_V_at)$(OCAMLOPT) $(OCAMLFLAGS) $(OCAMLOPTFLAGS) -a -o $@ $^

omega_models.cma: $(OMEGA_MODELS_CMO)
	@if $(AM_V_P); then :; else echo "  OCAMLC   " $@; fi
	$(AM_V_at)$(OCAMLC) $(OCAMLFLAGS) -a -o $@ $^

cascade_lexer.mli: cascade_lexer.ml cascade_parser.cmi
	$(OCAMLC) -i $< | $(GREP) 'val token' >$@

vertex_lexer.mli: vertex_lexer.ml vertex_parser.cmi
	$(OCAMLC) -i $< | $(EGREP) 'val (token|init_position)' >$@

UFO_lexer.mli: UFO_lexer.ml UFO_parser.cmi UFO_tools.cmi
	$(OCAMLC) -i $< | $(EGREP) 'val (token|init_position)' >$@

UFOx_lexer.mli: UFOx_lexer.ml UFOx_parser.cmi UFO_tools.cmi
	$(OCAMLC) -i $< | $(EGREP) 'val (token|init_position)' >$@

endif

MYPRECIOUS = $(OMEGA_DERIVED_CAML)

SUFFIXES += .lo .$(FCMOD)

# Fortran90 module files are generated at the same time as object files
.lo.$(FCMOD):
	@:
#	touch $@

########################################################################

DISTCLEANFILES = kinds.f90

if NOWEB_AVAILABLE

omegalib.stamp: $(srcdir)/omegalib.nw
	@rm -f omegalib.tmp
	@touch omegalib.tmp
	for src in $(OMEGALIB_DERIVED_F90); do \
	  $(NOTANGLE) -R[[$$src]] $< | $(CPIF) $$src; \
	done
	@mv -f omegalib.tmp omegalib.stamp

$(OMEGALIB_DERIVED_F90): omegalib.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
	  rm -f omegalib.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) omegalib.stamp; \
	fi

DISTCLEANFILES += $(OMEGALIB_DERIVED_F90)

endif NOWEB_AVAILABLE

MYPRECIOUS += $(OMEGALIB_DERIVED_F90)

########################################################################
# The following line just says
#    include Makefile.depend_fortran
# but in a portable fashion (depending on automake's AM_MAKE_INCLUDE
########################################################################
@am__include@ @am__quote@Makefile.depend_fortran@am__quote@

Makefile.depend_fortran: kinds.f90 $(libomega_core_la_SOURCES)
	@rm -f $@
	for src in $^; do \
	  module="`basename $$src | sed 's/\.f[90][0358]//'`"; \
	  grep '^ *use ' $$src \
	    | grep -v '!NODEP!' \
	    | sed -e 's/^ *use */'$$module'.lo: /' \
	          -e 's/, *only:.*//' \
	          -e 's/, *&//' \
	          -e 's/, *.*=>.*//' \
	          -e 's/ *$$/.lo/' ; \
	done > $@

DISTCLEANFILES += Makefile.depend_fortran 

if OCAML_AVAILABLE

@am__include@ @am__quote@Makefile.depend_ocaml@am__quote@

PARSERS = cascade vertex UFO UFOx

Makefile.depend_ocaml: $(OMEGA_CAML_PRIMARY)
	@if $(AM_V_P); then :; else echo "  OCAMLDEP " $@; fi
	@rm -f $@
	$(AM_V_at)$(OCAMLDEP) -I $(srcdir) $^ $(OMEGA_DERIVED_CAML) \
	  | sed 's,[^ 	]*/,,g' > $@
	$(AM_V_at)for parser in $(PARSERS); do \
	  echo $${parser}.cmi: $${parser}_lexer.cmi; \
	  echo $${parser}_lexer.cmi: $${parser}_parser.cmi; \
	  echo $${parser}_parser.cmi: $${parser}_syntax.cmi; \
	  echo $${parser}_parser.mli: $${parser}_parser.ml; \
	  echo $${parser}.cmo: $${parser}.cmi; \
	  echo $${parser}.cmx: $${parser}.cmi $${parser}_lexer.cmx; \
	  echo $${parser}_lexer.cmo: $${parser}_lexer.cmi; \
	  echo $${parser}_lexer.cmx: $${parser}_lexer.cmi $${parser}_parser.cmx; \
	  echo $${parser}_parser.cmo: $${parser}_parser.cmi $${parser}_syntax.cmi; \
	  echo $${parser}_parser.cmx: $${parser}_parser.cmi \
	    $${parser}_syntax.cmi $${parser}_syntax.cmx; \
	done >>$@

DISTCLEANFILES += Makefile.depend_ocaml

endif OCAML_AVAILABLE

########################################################################
# Don't trigger remakes by deleting intermediate files.
.PRECIOUS = $(MYPRECIOUS)

clean-local:
	rm -f *.cm[aiox] *.cmxa *.[ao] *.l[oa] *.$(FCMOD) \
		$(OMEGA_DERIVED_CAML) omegalib.stamp
if FC_SUBMODULES
	-rm -f *.smod
endif

distclean-local:
	-test "$(srcdir)" != "." && rm -f config.mli

########################################################################
## The End.
########################################################################
