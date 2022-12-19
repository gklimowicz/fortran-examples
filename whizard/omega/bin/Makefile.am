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

VPATH = $(top_builddir)/omega/src:$(top_srcdir)/omega/src

bin_SCRIPTS = $(OMEGA_APPLICATIONS_NATIVE)

bytecode: $(OMEGA_APPLICATIONS_BYTECODE)

OMEGA_CMXA = omega_core.cmxa omega_targets.cmxa omega_models.cmxa
OMEGA_CMA = $(OMEGA_CMXA:.cmxa=.cma)

########################################################################

include $(top_srcdir)/omega/src/Makefile.ocaml
include $(top_srcdir)/omega/src/Makefile.sources

OCAMLFLAGS += -I $(top_builddir)/omega/src

$(OMEGA_APPLICATIONS_NATIVE): $(OMEGA_CMXA)
$(OMEGA_APPLICATIONS_BYTECODE): $(OMEGA_CMA)

omegatop: $(OMEGA_CMA)
	$(OCAMLMKTOP) $(OCAMLFLAGS) -o $@ unix.cma $(OMEGA_CMA)

########################################################################

SUFFIXES += $(OCAML_NATIVE_EXT) .$(OMEGA_CACHE_SUFFIX) .$(OMEGA_CACHE_SUFFIX_COLORED)

$(OCAML_NATIVE_EXT).$(OMEGA_CACHE_SUFFIX):
	@echo $@:
	./$< -initialize .

.$(OMEGA_CACHE_SUFFIX).$(OMEGA_CACHE_SUFFIX_COLORED):
	@echo $@:

clean-local:
	rm -f *.cm[iox] *.o \
	  $(OMEGA_APPLICATIONS_BYTECODE) $(OMEGA_APPLICATIONS_NATIVE)

########################################################################
## The End.
########################################################################
