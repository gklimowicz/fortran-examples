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

TESTS = \
  circe2_tool.run \
  fail.run \
  flat.run moderate.run steep.run low.run \
  power.run resonance.run \
  multi.run delta.run \
  smooth.run oversmooth.run \
  scale.run 
# generate.run

XFAIL_TESTS =

BUILT_SOURCES = \
    OCAML_FLAG

OCAML_FLAG:
if OCAML_AVAILABLE
	touch $@
endif

EXTRA_DIST = $(TESTS:.run=.sh)

VPATH = $(srcdir)
SUFFIXES = .sh .run

.sh.run:
	@rm -f $@
	$(SED) 's|@name@|$*|g' $< >$@
	chmod +x $@

## Remove generated files
clean-local:
	rm -f *.output *.run *.circe2 OCAML_FLAG

maintainer-clean-local:
	-rm -f *~
.PHONY: maintainer-clean-local

if OCAML_AVAILABLE
installcheck-local:
	PATH=$(DESTDIR)$(bindir):$$PATH; export PATH; \
	LD_LIBRARY_PATH=$(DESTDIR)$(libdir):$(DESTDIR)$(pkglibdir):$$LD_LIBRARY_PATH; \
	export LD_LIBRARY_PATH; \
	circe2_tool.opt -test || circe2_tool.bin -test || circe2_tool -test
endif
