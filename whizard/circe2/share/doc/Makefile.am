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

DISTCLEANFILES =

include $(top_srcdir)/circe2/src/Makefile.sources

VPATH = $(srcdir):$(top_builddir)/circe2/src:$(srcdir):$(top_srcdir)/circe2/src

WEBS = prelude.nw circe2.nw postlude.nw

if DISTRIBUTION
PDFS = circe2.pdf 
else
PDFS = 
endif

LATEX_STYLES = \
    emp.sty noweb.sty thohacks.sty thopp.sty ocamlweb.sty

TEX_FLAGS = "$$TEXINPUTS:$(top_srcdir)/circe2/share/doc"
MP_FLAGS  = "$$MPINPUTS:$(top_srcdir)/circe2/share/doc"

CIRCE2_HISTOS = \
  x.20.histo x.20m.histo x.20q.histo x.20qm.histo x.input.histo \
  z.20.histo z.20m2.histo z.20m.histo z.20q.histo z.20qm.histo \
  z.50m2.histo z.input.histo \
  z_low.20.histo z_low.20m2.histo z_low.20m.histo z_low.20q.histo \
  z_low.20qm.histo z_low.50.histo z_low.50m2.histo z_low.50q.histo \
  z_low.input.histo

CIRCE2_HISTOSDATA = $(CIRCE2_HISTOS:.histo=.data)

DISTCLEANFILES += $(CIRCE2_HISTOS)

EXTRA_DIST = \
  tex-comments.sh \
  $(LATEX_STYLES) \
  $(CIRCE2_HISTOSDATA)

if DISTRIBUTION
dist_doc_DATA = $(PDFS)
endif

if NOWEB_AVAILABLE
pdf-local: circe2.pdf 
endif

if NOWEB_AVAILABLE

circe2.tex: $(WEBS)
	@if $(AM_V_P); then :; else echo " NOWEAVE   " $@; fi
	$(AM_V_at)$(NOWEAVE) -filter ./tex-comments -delay -index \
	    `for i in  $^; do case $$i in *.nw) echo $$i;; esac done` \
	  > $@

circe2.tex: tex-comments

endif NOWEB_AVAILABLE

tex-comments: tex-comments.sh
	cp $< $@
	chmod +x $@

.data.histo:
	@if $(AM_V_P); then :; else echo " CP        " $@; fi
	$(AM_V_at)cp $< $@

# preview.pdf: vegas.data vamp.data

# vegas.data: vegas.d
# 	cp $< $@
# 
# vamp.data: vamp.d
# 	cp $< $@

SUFFIXES = \
  .mly .mll .ml .implementation .mli .interface \
  .data .histo .tex .pdf \
  .nw .dvi .eps .ps

if !NOWEB_AVAILABLE

circe2.pdf:

else NOWEB_AVAILABLE

circe2.pdf: $(CIRCE2_INTERFACE) $(CIRCE2_IMPLEMENTATION) \
		$(CIRCE2TOOL_IMPLEMENTATION) circe2.tex \
		index.tex $(CIRCE2_HISTOS)

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

index.tex: $(CIRCE2_CAML) $(CIRCE2_DERIVED)
	@if $(AM_V_P); then :; else echo " OCAMLWEB  " $@; fi
	$(AM_V_at)$(OCAMLWEB) --no-preamble --noweb $^ | \
	   sed -n '/\\ocwbeginindex{}/,/\\ocwendindex{}/p' >$@

endif OCAMLWEB_AVAILABLE
endif DISTRIBUTION

endif NOWEB_AVAILABLE

MPOST_LATEX = TEX=$(LATEX) TEXINPUTS=$(TEX_FLAGS) $(MPOST)

if DISTRIBUTION
if PDFLATEX_AVAILABLE
if SUPP_PDF_AVAILABLE
.tex.pdf:
	@if $(AM_V_P); then TEXINPUTS=$(TEX_FLAGS) $(PDFLATEX) $<; else \
	   echo " PDFLATEX  " $<; TEXINPUTS=$(TEX_FLAGS) $(PDFLATEX) $< >/dev/null; fi
if MAKEINDEX_AVAILABLE
	@if $(AM_V_P); then $(MAKEINDEX) -o $*.ind $*.idx; else \
	   echo " MAKEINDEX " $*.ind $*.idx; $(MAKEINDEX) -q -o $*.ind $*.idx; fi
endif MAKEINDEX_AVAILABLE
if MPOST_AVAILABLE
	@if $(AM_V_P); then test -r $*.mp && $(MPOST_LATEX) $*; else \
	   echo " METAPOST  " $*.mp; test -r $*.mp && $(MPOST_LATEX) $* >/dev/null; fi
endif MPOST_AVAILABLE
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
endif SUPP_PDF_AVAILABLE
endif PDFLATEX_AVAILABLE
endif DISTRIBUTION

## Cleanup tasks
mostlyclean-latex:
	-rm -f *.mpx *.[1-9]* *.t[1-9]* circe*.mp preview*.mp \
	   *.out *.log *.aux *.idx *.ilg *.ind tex-comments *.toc \
	   circe2.tex
clean-latex:
maintainer-clean-latex:
	-rm $(PDFS)
if NOWEB_AVAILABLE
if OCAMLWEB_AVAILABLE
mostlyclean-circe2:
	-test "$(srcdir)" != "." && rm -f $(PDFS)
maintainer-clean-circe2:
else
mostlyclean-circe2:
maintainer-clean-circe2:
endif
endif
.PHONY: mostlyclean-latex clean-latex maintainer-clean-latex
.PHONY: mostlyclean-circe2 maintainer-clean-circe2

if OCAMLWEB_AVAILABLE
mostlyclean-caml:
	-rm -f *.interface *.implementation index.tex
else
mostlyclean-caml:
endif
clean-caml:
if OCAMLWEB_AVAILABLE
maintainer-clean-caml:
	-rm -f *.interface *.implementation index.tex
else
maintainer-clean-caml:
endif
.PHONY: mostlyclean-caml clean-caml maintainer-clean-caml

## Remove backup files
maintainer-clean-backup:
	-rm -f *~
.PHONY: maintainer-clean-backup

## Register additional clean targets
mostlyclean-local: mostlyclean-latex mostlyclean-circe2 \
		mostlyclean-caml
clean-local: clean-latex clean-caml
maintainer-clean-local: maintainer-clean-latex maintainer-clean-circe2 \
		maintainer-clean-caml maintainer-clean-backup

if !DISTRIBUTION
install-data-hook: 
	-$(INSTALL) -m 644 circe2.pdf $(DESTDIR)$(datarootdir)/doc/circe2

uninstall-hook: 
	-rm -f $(DESTDIR)/$(datarootdir)/doc/circe2/circe2.pdf 
endif

########################################################################
## The End.
########################################################################

