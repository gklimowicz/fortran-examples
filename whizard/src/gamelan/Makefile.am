## Makefile.am -- Makefile for WHIZARD
##
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
## Process this file with automake to produce Makefile.in

## Install the gml executable and WHIZARD logo here:
GML=whizard-gml
gmlexeclibdir = $(pkglibdir)/gamelan
gmlexecbindir = $(bindir)
dist_gmlexecbin_SCRIPTS = $(GML)

## Install the gml MetaPost code here:
gmllibdir = $(pkglibdir)/gamelan

## The NOWEB source and the test / demofiles generated from it must be distributed,
## but not installed.
dist_noinst_DATA = \
  gamelan.nw \
  $(GML_TEX_EXAMPLES_SRC) $(GML_TEX_EXAMPLES_DAT) \
  $(GML_TEX_TEST_SRC) $(GML_TEX_TEST_DAT)

## Mem files are no longer existent in MetaPost, hence we need to install 
## all the gamelan includes.
dist_gmllib_DATA = $(GML_MPOST_SRC) $(GML_WHIZARD_LOGO)

## Distribute the style file and color definition files, and install them here:
gmltexdir = $(PKGTEXDIR)
dist_gmltex_DATA = \
  $(GML_STY) \
  $(GML_TEX_SRC)

## ---------------------------------------------------------------
## File lists

## Contents of gamelan.nw
GML_FROM_GAMELAN_NW = \
  $(GML_MPOST_SRC) $(GML_STY) $(GML_TEX_SRC) \
  $(GML_TEX_EXAMPLES_SRC) $(GML_TEX_EXAMPLES_DAT) \
  $(GML_TEX_TEST_SRC) $(GML_TEX_TEST_DAT) $(GML_WHIZARD_LOGO)

## The WHIZARD logo as EPS file
GML_WHIZARD_LOGO = \
  Whizard-Logo.eps

## Gamelan MetaPost sources
GML_MPOST_SRC = \
  gamelan.mp \
  gmlaux.mp gmlhatch.mp gmlbox.mp gmlshapes.mp gmlspectra.mp \
  gmlarith.mp gmlformat.mp gmlfiles.mp gmlset.mp gmlscan.mp \
  gmlcoords.mp gmldraw.mp gmlgrid.mp gmleps.mp

## LaTeX style file
GML_STY = \
  gamelan.sty

## Color definitions, to be included upon request:
GML_TEX_SRC = \
  gmlcolors.tex gmlextracolors.tex

## Examples for gamelan usage:
GML_TEX_EXAMPLES_SRC = \
  gmlfun.tex gmldata.tex gmlerr.tex gmlimp.tex
GML_TEX_EXAMPLES_DAT = \
  gmldata.dat gmlerr.dat 
GML_TEX_EXAMPLES_MP = \
  gmlfun.mp gmldata.mp gmlerr.mp gmlimp.mp
GML_TEX_EXAMPLES_TARGETS = \
  gmlfun.ps gmldata.ps gmlerr.ps gmlimp.ps

## Self-tests
GML_TEX_TEST_SRC = \
  gmltest.tex gmllongtest.tex 
GML_TEX_TEST_DAT = \
  gmllongtest.dat 
GML_TEX_TEST_MP = \
  gmltest.mp gmllongtest.mp
GML_TEX_TEST_TARGETS = 
if DVIPS_AVAILABLE
GML_TEX_TEST_TARGETS += \
  gmltest.ps gmllongtest.ps
endif

## --------------------------------------------------------------------
## Special targets

## The 'check' rule consists of building these PS files from TEX sources:
check: $(GML_TEX_TEST_TARGETS)

## These examples can be built upon request:
examples: $(GML_TEX_EXAMPLES_TARGETS)

## don't try to run the examples in parallel (TeXLive 2009 doesn't like it)
gmllongtest.ps: gmltest.ps

## Run texhash to include the style file in the TeX search path
install-data-hook:
	test -x `which texhash` && texhash

## Create sources from the noweb file, if possible
if NOWEB_AVAILABLE
gamelan.stamp: $(srcdir)/gamelan.nw
	@rm -f gamelan.tmp
	@touch gamelan.tmp
	for src in $(GML_FROM_GAMELAN_NW); do \
	  $(NOTANGLE) -R$$src $< | $(CPIF) $$src; \
        done
	@mv -f gamelan.tmp gamelan.stamp

$(GML_FROM_GAMELAN_NW): gamelan.stamp
## Recover from the removal of $@
	@if test -f $@; then :; else \
        rm -f gamelan.stamp; \
	  $(MAKE) $(AM_MAKEFLAGS) gamelan.stamp; \
	fi

endif

## Implicit rule for building test and examples
## If noweb is not available copy the stuff over
if NOWEB_AVAILABLE
.tex.ps: 
	$(LATEX) $*
	./whizard-gml $*
	$(LATEX) $*
	$(DVIPS) -o $@ $*
else !NOWEB_AVAILABLE
.tex.ps:
	-if test "$(srcdir)" != "."; then \
          for file in $(GML_FROM_GAMELAN_NW); do \
            test -f "$$file" || cp $(srcdir)/$$file .; \
          done; \
        fi
	$(LATEX) $*
	./whizard-gml $*
	$(LATEX) $*
	$(DVIPS) -o $@ $*
endif

## Non-standard cleanup tasks
if NOWEB_AVAILABLE
maintainer-clean-noweb:
	-rm -f $(GML_FROM_GAMELAN_NW)
endif
.PHONY: maintainer-clean-noweb

## Remove these sources also if builddir and srcdir are different
if NOWEB_AVAILABLE
clean-noweb:
	-rm -f gamelan.stamp gamelan.tmp
	test "$(srcdir)" != "." && rm -f $(GML_FROM_GAMELAN_NW) || true
	-rm -f $(GML_TEX_EXAMPLES_MP) $(GML_TEX_TEST_MP)
else !NOWEB_AVAILABLE
clean-noweb:
	test "$(srcdir)" != "." && rm -f *.sty *.tex *.mp *.dat *.eps || true
endif
.PHONY: clean-noweb

## Remove test 

clean-local: clean-noweb
	-rm -f *.ps *.pool *.aux *.dvi *.idx *.log *.toc
	-rm -f *.mpx *.sh *.tmp *.1 *.ltp 

## Remove backup files
maintainer-clean-local: maintainer-clean-noweb
	-rm -f *~ 
