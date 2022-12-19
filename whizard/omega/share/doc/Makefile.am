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

########################################################################
### TODO: fix weaving of lexers and parsers
########################################################################

include $(top_srcdir)/omega/src/Makefile.sources

VPATH = $(srcdir):$(top_builddir)/omega/src:$(srcdir):$(top_srcdir)/omega/src

PICTURES_PDF = \
    modules.pdf \
    omega-paper-1-pics-1.pdf \
    omega-paper-1-pics-2.pdf \
    omega-paper-1-pics-3.pdf \
    omega-paper-1-pics-4.pdf \
    omega-paper-1-pics-5.pdf \
    omega-paper-1-pics-6.pdf \
    omega-paper-1-pics-7.pdf \
    omega-paper-1-pics-8.pdf \
    omega-paper-1-pics-9.pdf \
    omega-paper-1-pics-10.pdf \
    bhabha.pdf bhabha0.pdf \
    epemudbardubar.pdf epemudbardubar0.pdf \
    epemudbarmunumubar.pdf epemudbarmunumubar0.pdf \
    sign_ex.pdf fusion_rules.pdf mom_choice.pdf \
    mom_flow.pdf

LATEX_STYLES = \
    flex.cls thophys.sty thohacks.sty \
    noweb.sty ocamlweb.sty ytableau.sty \
    feynmp.sty feynmp.mp emp.sty

TEX_FLAGS = "$$TEXINPUTS:$(top_srcdir)/omega/share/doc"
MP_FLAGS  = "$$MPINPUTS:$(top_srcdir)/omega/share/doc"

if DISTRIBUTION
PDFS = omega.pdf omega-paper-1.pdf omega-paper-2.pdf 
else
PDFS = 
endif

### Files needed to be installed with the O'Mega distribution
modelsdir = $(pkgdatadir)/doc
if SUPP_PDF_AVAILABLE
dist_doc_DATA = $(PDFS)
else
dist_doc_DATA = 
endif

EXTRA_DIST = $(PICTURES_PDF) $(LATEX_STYLES)

if NOWEB_AVAILABLE
pdf-local: $(PDFS)
else
pdf-local:
endif

SUFFIXES = .mly .mll .ml .implementation .mli .interface .nw .tex .pdf

MPOST_LATEX = TEX=$(LATEX) TEXINPUTS=$(TEX_FLAGS) MPINPUTS=$(MP_FLAGS) $(MPOST)

if DISTRIBUTION
if SUPP_PDF_AVAILABLE
if PDFLATEX_AVAILABLE
.tex.pdf:
	@if $(AM_V_P); then TEXINPUTS=$(TEX_FLAGS) $(PDFLATEX) $<; else \
	   echo " PDFLATEX  " $<; TEXINPUTS=$(TEX_FLAGS) $(PDFLATEX) $< >/dev/null; fi
if MPOST_AVAILABLE
	@if test -r $*pics.mp; then \
	if $(AM_V_P); then MPINPUTS=$(MP_FLAGS) $(MPOST) $*pics; else \
	   echo " METAPOST  " $*pics.mp; MPINPUTS=$(MP_FLAGS) $(MPOST) $*pics > /dev/null; fi; \
	fi
	@if test -r $*.mp; then \
	if $(AM_V_P); then $(MPOST_LATEX) $*; else \
	   echo " METAPOST  " $*.mp; $(MPOST_LATEX) $* >/dev/null; fi; \
	fi
endif MPOST_AVAILABLE
	$(AM_V_at)echo " PDFLATEX   skipping -bibtex $*"
	@if $(AM_V_P); then TEXINPUTS=$(TEX_FLAGS) $(PDFLATEX) $<; else \
	   echo " PDFLATEX  " $<; TEXINPUTS=$(TEX_FLAGS) $(PDFLATEX) $< >/dev/null; fi
	@if $(AM_V_P); then \
	   if grep -s 'Rerun to get cross-references right.' $*.log; then \
	      TEXINPUTS=$(TEX_FLAGS) $(PDFLATEX) $<; \
	   fi; else \
	   if grep -s 'Rerun to get cross-references right.' $*.log >/dev/null; then \
	      echo " PDFLATEX  " $< "(for cross-references)"; \
	      TEXINPUTS=$(TEX_FLAGS) $(PDFLATEX) $< >/dev/null; \
	   fi; \
	fi
endif PDFLATEX_AVAILABLE
endif SUPP_PDF_AVAILABLE
endif DISTRIBUTION

if DISTRIBUTION
if SUPP_PDF_AVAILABLE
if PDFLATEX_AVAILABLE

omega-paper-1.pdf: modules.pdf \
    omega-paper-1-pics-1.pdf \
    omega-paper-1-pics-2.pdf \
    omega-paper-1-pics-3.pdf \
    omega-paper-1-pics-4.pdf \
    omega-paper-1-pics-5.pdf \
    omega-paper-1-pics-6.pdf \
    omega-paper-1-pics-7.pdf \
    omega-paper-1-pics-8.pdf \
    omega-paper-1-pics-9.pdf \
    omega-paper-1-pics-10.pdf

# Dependencies and avoid mpost race condition
omega-paper-2.pdf: \
    omega-paper-1.pdf sign_ex.pdf fusion_rules.pdf \
    mom_choice.pdf mom_flow.pdf

endif PDFLATEX_AVAILABLE
endif SUPP_PDF_AVAILABLE
endif DISTRIBUTION

OMEGA_CORE_INTERFACES = $(OMEGA_CORE_MLI:.mli=.interface)
OMEGA_CORE_IMPLEMENTATIONS = $(OMEGA_CORE_ML:.ml=.implementation)
OMEGA_MODELLIB_INTERFACES = $(OMEGA_MODELLIB_MLI:.mli=.interface)
OMEGA_MODELLIB_IMPLEMENTATIONS = $(OMEGA_MODELLIB_ML:.ml=.implementation)
OMEGA_TARGETLIB_INTERFACES = $(OMEGA_TARGETLIB_MLI:.mli=.interface)
OMEGA_TARGETLIB_IMPLEMENTATIONS = $(OMEGA_TARGETLIB_ML:.ml=.implementation)
OMEGA_APPLICATIONS_IMPLEMENTATIONS = $(OMEGA_APPLICATIONS_ML:.ml=.implementation)

OMEGA_INTERFACES = \
    $(OMEGA_CORE_INTERFACES) \
    $(OMEGA_MODELLIB_INTERFACES) \
    $(OMEGA_TARGETLIB_INTERFACES)

OMEGA_IMPLEMENTATIONS = \
    $(OMEGA_CORE_IMPLEMENTATIONS) \
    $(OMEGA_MODELLIB_IMPLEMENTATIONS) \
    $(OMEGA_TARGETLIB_IMPLEMENTATIONS) \
    $(OMEGA_APPLICATIONS_IMPLEMENTATIONS)

if !NOWEB_AVAILABLE

omega.pdf:

else NOWEB_AVAILABLE

omega.pdf: \
    $(OMEGA_INTERFACES) $(OMEGA_IMPLEMENTATIONS) omegalib.tex index.tex \
    $(PICTURES_PDF)

.nw.tex:
	@if $(AM_V_P); then :; else echo " NOWEAVE   " $@; fi
	$(AM_V_at)$(NOWEAVE) -delay $< > $@

if DISTRIBUTION
if OCAMLWEB_AVAILABLE

.mll.implementation:
	@if $(AM_V_P); then :; else echo " OCAMLWEB  " $@; fi
	$(AM_V_at)$(OCAMLWEB) --no-preamble --noweb --no-index $< >$@

.mly.implementation:
	@if $(AM_V_P); then :; else echo " OCAMLWEB  " $@; fi
	$(AM_V_at)$(OCAMLWEB) --no-preamble --noweb --no-index $< >$@

.ml.implementation:
	@if $(AM_V_P); then :; else echo " OCAMLWEB  " $@; fi
	$(AM_V_at)$(OCAMLWEB) --no-preamble --noweb --no-index $< >$@

.mli.interface:
	@if $(AM_V_P); then :; else echo " OCAMLWEB  " $@; fi
	$(AM_V_at)$(OCAMLWEB) --no-preamble --noweb --no-index $< >$@

index.tex: $(OMEGA_CAML)
	@if $(AM_V_P); then :; else echo " OCAMLWEB  " $@; fi
	$(AM_V_at)$(OCAMLWEB) --no-preamble --noweb $^ | \
	sed -n '/\\ocwbeginindex{}/,/\\ocwendindex{}/p' >$@

endif OCAMLWEB_AVAILABLE
endif DISTRIBUTION

endif NOWEB_AVAILABLE

########################################################################

## Cleanup tasks
mostlyclean-latex:
	-rm -f *.log *.aux *.toc *.mpx *.idx *.out omega*.mp \
	omega*pics.t[0-9]* omega*pics.[0-9]* $(PICTURES_PDF) \
	omegalib.tex
clean-latex:
maintainer-clean-latex:
	-rm $(PDFS)
if NOWEB_AVAILABLE
mostlyclean-omega:	
	-test "$(srcdir)" != "." && rm -f $(PDFS)
maintainer-clean-omega:
else
mostlyclean-omega:
maintainer-clean-omega:
endif 
.PHONY: mostlyclean-latex clean-latex maintainer-clean-latex
.PHONY: mostlyclean-omega maintainer-clean-omega

if OCAMLWEB_AVAILABLE
mostlyclean-caml:
	-rm -f $(OMEGA_INTERFACES) $(OMEGA_IMPLEMENTATIONS) index.tex
else
mostlyclean-caml:
endif
clean-caml:
if OCAMLWEB_AVAILABLE
maintainer-clean-caml:
	-rm -f $(OMEGA_INTERFACES) $(OMEGA_IMPLEMENTATIONS) index.tex
else 
maintainer-clean-caml:
endif
.PHONY: mostlyclean-caml clean-caml maintainer-clean-caml

## Remove backup files
maintainer-clean-backup:
	-rm -f *~
.PHONY: maintainer-clean-backup

## Register additional clean targets
mostlyclean-local: mostlyclean-latex mostlyclean-caml mostlyclean-omega
clean-local: clean-latex clean-caml
maintainer-clean-local: maintainer-clean-latex maintainer-clean-caml \
	maintainer-clean-omega maintainer-clean-backup

if !DISTRIBUTION
install-data-hook: 
	-$(INSTALL) -m 644 omega.pdf $(DESTDIR)$(datarootdir)/doc/omega
	-$(INSTALL) -m 644 omega-paper-1.pdf $(DESTDIR)$(datarootdir)/doc/omega
	-$(INSTALL) -m 644 omega-paper-2.pdf $(DESTDIR)$(datarootdir)/doc/omega

uninstall-hook: 
	-rm -f $(DESTDIR)/$(datarootdir)/doc/omega/omega.pdf 
	-rm -f $(DESTDIR)/$(datarootdir)/doc/omega/omega-paper-1.pdf 
	-rm -f $(DESTDIR)/$(datarootdir)/doc/omega/omega-paper-2.pdf 
endif

########################################################################
## The End.
########################################################################
