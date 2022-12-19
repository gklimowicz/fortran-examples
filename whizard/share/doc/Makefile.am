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

## The WHIZARD documented source is assembled from various directories
## defined outside the DISTRIBUTION environment as modern autotools
## versions complain otherwise
w_srcdir = $(top_srcdir)/src
VPATH = $(srcdir):$(w_srcdir)/noweb-frame:$(w_srcdir)/utilities:$(w_srcdir)/testing:$(w_srcdir)/system:$(w_srcdir)/combinatorics:$(w_srcdir)/parsing:$(w_srcdir)/rng:$(w_srcdir)/expr_base:$(w_srcdir)/physics:$(w_srcdir)/qed_pdf:$(w_srcdir)/qft:$(w_srcdir)/types:$(w_srcdir)/matrix_elements:$(w_srcdir)/particles:$(w_srcdir)/beams:$(w_srcdir)/me_methods:$(w_srcdir)/events:$(w_srcdir)/phase_space:$(w_srcdir)/vegas:$(w_srcdir)/mci:$(w_srcdir)/fks:$(w_srcdir)/gosam:$(w_srcdir)/openloops:$(w_srcdir)/blha:$(w_srcdir)/shower:$(w_srcdir)/muli:$(w_srcdir)/variables:$(w_srcdir)/model_features:$(w_srcdir)/threshold:$(w_srcdir)/process_integration:$(w_srcdir)/matching:$(w_srcdir)/transforms:$(w_srcdir)/whizard-core:$(w_srcdir)/main:$(w_srcdir)/api

## The primary targets
if DISTRIBUTION
## The manual source has to be distributed
dist_noinst_DATA = manual.tex $(PACKAGE).tex \
	book.hva custom.hva fancysection.hva Whizard-Logo.jpg \
	$(MANUAL_PICS) dep2dot.py

MANUAL_PICS = \
   proc_4f-history.pdf whizstruct.pdf cc10_1.pdf \
   cc10_2.pdf Z-lineshape_1.pdf Z-lineshape_2.pdf \
   flow4.pdf lep_higgs_1.pdf \
   lep_higgs_2.pdf lep_higgs_3.pdf circe2-smoothing.pdf \
   resonance_e_gam.pdf resonance_n_charged.pdf \
   resonance_n_hadron.pdf resonance_n_particles.pdf  \
   resonance_n_photons.pdf resonance_n_visible.pdf

if NOWEB_AVAILABLE
dist_pdf_DATA = manual.pdf $(PACKAGE).pdf gamelan_manual.pdf
else
dist_pdf_DATA = manual.pdf gamelan_manual.pdf
endif
else
dist_pdf_DATA =
endif

pdf-local: manual.pdf $(PACKAGE).pdf gamelan_manual.pdf

if DISTRIBUTION
if HEVEA_AVAILABLE
html_DATA = manual.html index.html
endif HEVEA_AVAILABLE
endif DISTRIBUTION

GML=../../src/gamelan/whizard-gml --math=scaled -halt-on-error --gmldir ../../src/gamelan

LATEX_STYLES = \
    noweb.sty thophys.sty gamelan.sty hevea.sty

TEX_FLAGS = "$$TEXINPUTS:$(top_srcdir)/share/doc"

EXTRA_DIST = $(LATEX_STYLES)

## don't try to run the files in parallel (TeXLive 2009 doesn't like it)
manual.pdf: $(PACKAGE).pdf
$(PACKAGE).pdf: gamelan_manual.pdf $(WHIZARD_DEPENDENCY_GRAPHS_PDF) overview.pdf

manual.pdf: variables.tex

variables.tex: ../../src/whizard
	../../src/whizard --generate-variables-tex > variables.tex

gamelan_manual.pdf: gamelan.sty

SUFFIXES: .dot .tex .pdf

.dot.pdf:
	@if $(AM_V_P); then :; else echo " DOT       " $@; fi
	$(AM_V_at)$(DOT) -Tpdf $< > $@

if DISTRIBUTION
WHIZARD_DEPENDENCY_GRAPHS_DOT = \
	blha.dot \
	beams.dot \
	combinatorics.dot \
	events.dot \
	expr_base.dot \
	fks.dot \
	gosam.dot \
	matching.dot \
	matrix_elements.dot \
	vegas.dot \
	mci.dot \
	me_methods.dot \
	model_features.dot \
	muli.dot \
	openloops.dot \
	parsing.dot \
	particles.dot \
	phase_space.dot \
	physics.dot \
	qed_pdf.dot \
	process_integration.dot \
	qft.dot \
	rng.dot \
	shower.dot \
	system.dot \
	testing.dot \
	threshold.dot \
	transforms.dot \
	types.dot \
	utilities.dot \
	variables.dot \
	whizard-core.dot \
	main.dot \
	api.dot

WHIZARD_DEPENDENCY_GRAPHS_PDF = $(WHIZARD_DEPENDENCY_GRAPHS_DOT:.dot=.pdf)
all-dots: $(WHIZARD_DEPENDENCY_GRAPHS_DOT) overview.dot

$(WHIZARD_DEPENDENCY_GRAPHS_DOT):
	@rm -f $@
	@if $(AM_V_P); then \
	   $(top_srcdir)/share/doc/dep2dot.py \
		../../src/`echo $@ | sed 's/.dot//'`/Makefile.depend > $@; else \
	   echo " DEP2DOT   " $@; \
	   $(top_srcdir)/share/doc/dep2dot.py \
		../../src/`echo $@ | sed 's/.dot//'`/Makefile.depend > $@; fi

overview.dot:
	@rm -f $@
	@if $(AM_V_P); then \
	   list=''; \
	   for dep in $(WHIZARD_DEPENDENCY_GRAPHS_DOT); do \
		list="$$list ../../src/`echo $$dep | sed 's/.dot//'`/Makefile.depend"; \
	   done ; \
	   $(top_srcdir)/share/doc/dep2dot.py $$list > $@; else \
	   echo " DEP2DOT   " $@; \
	   list=''; \
	   for dep in $(WHIZARD_DEPENDENCY_GRAPHS_DOT); do \
		list="$$list ../../src/`echo $$dep | sed 's/.dot//'`/Makefile.depend"; \
	   done ; \
	   $(top_srcdir)/share/doc/dep2dot.py $$list > $@; fi
dist_noinst_DATA += $(WHIZARD_DEPENDENCY_GRAPHS_DOT) overview.dot
endif

## Rules for creating PDF
if DISTRIBUTION
if PDFLATEX_AVAILABLE
.tex.pdf:
	@if $(AM_V_P); then TEXINPUTS=$(TEX_FLAGS) $(PDFLATEX) $<; else \
	   echo " PDFLATEX  " $<; TEXINPUTS=$(TEX_FLAGS) $(PDFLATEX) $< >/dev/null; fi
	@if $(AM_V_P); then \
	   while grep 'Rerun to get cross-references right\.' $*.log; \
	   do TEXINPUTS=$(TEX_FLAGS) $(PDFLATEX) $<; \
	   done; else \
	   echo " PDFLATEX  " $< "(for cross-references)"; \
	   while grep 'Rerun to get cross-references right\.' $*.log >/dev/null; \
	   do TEXINPUTS=$(TEX_FLAGS) $(PDFLATEX) $< >/dev/null; \
	   done; \
	fi
	@if test -r $*.mp; then \
	   if $(AM_V_P); then $(GML) $*; else echo " GML       " $*; $(GML) $* >/dev/null; fi; \
	fi
	@if $(AM_V_P); then TEXINPUTS=$(TEX_FLAGS) $(PDFLATEX) $<; else \
	   echo " PDFLATEX  " $<; TEXINPUTS=$(TEX_FLAGS) $(PDFLATEX) $< >/dev/null; fi
endif PDFLATEX_AVAILABLE
endif DISTRIBUTION

## Rules for creating HTML
if HEVEA_AVAILABLE
HEVEAOPTS = -exec xxdate.exe -I $(top_srcdir)/share/doc \
	book.hva fancysection.hva custom.hva
HACHAOPTS = -tocbis

if DISTRIBUTION
index.html: manual.html
	@if $(AM_V_P); then $(HACHA) $(HACHAOPTS) -o index.html manual.html; else \
	   echo " HACHA     " $@; $(HACHA) $(HACHAOPTS) -o index.html manual.html >/dev/null 2>&1; fi

manual.html: variables.tex manual.tex $(MANUAL_PICS)
	@if $(AM_V_P); then TEXINPUTS=$(TEX_FLAGS) $(HEVEA) -fix $(HEVEAOPTS) manual.tex; else \
	   echo " HEVEA     " $@; \
	   TEXINPUTS=$(TEX_FLAGS) $(HEVEA) -s -fix $(HEVEAOPTS) manual.tex >/dev/null 2>&1; fi
	@if $(AM_V_P); then TEXINPUTS=$(TEX_FLAGS) $(IMAGEN) -pdf manual; else \
	   echo " IMAGEN     manual"; TEXINPUTS=$(TEX_FLAGS) $(IMAGEN) -pdf manual >/dev/null 2>&1; fi
### There are no Feynman diagrams at the moment inside the manual.
#	$(MPOST) manualpics.mp
	@if $(AM_V_P); then TEXINPUTS=$(TEX_FLAGS) $(HEVEA) -fix $(HEVEAOPTS) manual.tex; else \
	   echo " HEVEA     " $@; TEXINPUTS=$(TEX_FLAGS) $(HEVEA) -s -fix $(HEVEAOPTS) manual.tex; fi
else
	@echo "HEVEA not available. The HTML manual cannot be made"
endif !DISTRIBUTION
endif HEVEA_AVAILABLE

if NOWEB_AVAILABLE
WHIZARD_NOWEB_SRC = \
  whizard-prelude.nw \
  utilities.nw \
  testing.nw \
  system.nw \
  combinatorics.nw \
  parsing.nw rng.nw physics.nw qed_pdf.nw qft.nw \
  types.nw \
  matrix_elements.nw \
  particles.nw \
  beams.nw \
  me_methods.nw \
  events.nw \
  phase_space.nw \
  vegas.nw \
  mci.nw \
  shower.nw muli.nw \
  blha.nw gosam.nw \
  openloops.nw fks.nw \
  model_features.nw \
  threshold.nw \
  process_integration.nw \
  matching.nw \
  transforms.nw \
  whizard.nw \
  whizard-postlude.nw

$(PACKAGE).tex: $(WHIZARD_NOWEB_SRC)
	-rm -f $@
	@if $(AM_V_P); then :; else echo " NOWEAVE   " $@; fi
	$(AM_V_at)$(NOWEAVE) -delay $^ | $(CPIF) $@

gamelan.sty: $(top_srcdir)/src/gamelan/gamelan.nw
	@if $(AM_V_P); then :; else echo " NOTANGLE  " $@; fi
	$(AM_V_at)$(NOTANGLE) -R$@ $< | $(CPIF) $@
endif NOWEB_AVAILABLE

## Cleanup tasks
mostlyclean-latex:
	-rm -f *.aux *.log *.dvi *.toc *.idx *.out *.ltp *.mp *.mpx *.glo \
	gamelan_manual.[1-9] gamelan_manual.[1-9][0-9] \
	manual.pdf gamelan_manual.pdf \
	$(WHIZARD_DEPENDENCY_GRAPHS_DOT) $(WHIZARD_DEPENDENCY_GRAPHS_PDF) \
	overview.dot overview.pdf variables.tex
	-test "$(srcdir)" != "." && rm -f $(PACKAGE).pdf
clean-latex:
maintainer-clean-latex:
	-rm manual.pdf gamelan_manual.pdf
if NOWEB_AVAILABLE
mostlyclean-whizard:
	-rm -f $(PACKAGE).tex $(PACKAGE).pdf gamelan.sty
maintainer-clean-whizard:
else
mostlyclean-whizard:
maintainer-clean-whizard:
endif
.PHONY: mostlyclean-latex clean-latex maintainer-clean-latex
.PHONY: mostlyclean-whizard maintainer-clean-whizard

mostlyclean-html:
	-rm -f *.haux *.htoc *.css index.html contents_motif.gif \
	next_motif.gif previous_motif.gif contents_motif.svg \
	next_motif.svg previous_motif.svg manual*.html  manual*.png \
	manual.image.tex
clean-html:
maintainer-clean-html:
	-rm -f manual*.html index.html contents_motif.gif \
	next_motif.gif previous_motif.gif contents_motif.svg \
	next_motif.svg previous_motif.svg
.PHONY: mostlyclean-html clean-html maintainer-clean-html

## Remove backup files
maintainer-clean-backup:
	-rm -f *~
.PHONY: maintainer-clean-backup

## Register additional clean targets
mostlyclean-local: mostlyclean-latex mostlyclean-html mostlyclean-whizard
clean-local: clean-latex clean-html
maintainer-clean-local: maintainer-clean-latex maintainer-clean-html \
	maintainer-clean-whizard maintainer-clean-backup

if !DISTRIBUTION
install-data-hook:
	-$(INSTALL) -m 644 manual.pdf $(DESTDIR)$(datarootdir)/doc/whizard
	-$(INSTALL) -m 644 $(PACKAGE).pdf $(DESTDIR)$(datarootdir)/doc/whizard
	-$(INSTALL) -m 644 gamelan_manual.pdf $(DESTDIR)$(datarootdir)/doc/whizard

uninstall-hook:
	-rm -f $(DESTDIR)/$(datarootdir)/doc/whizard/manual.pdf
	-rm -f $(DESTDIR)/$(datarootdir)/doc/whizard/$(PACKAGE).pdf
	-rm -f $(DESTDIR)/$(datarootdir)/doc/whizard/gamelan_manual.pdf
endif

########################################################################
## The End.
########################################################################
