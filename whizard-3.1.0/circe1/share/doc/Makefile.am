# Makefile.am --
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

CIRCE1_VER = 1
CIRCE1_REV = 0

WEBS = \
  prelude.nw \
  circe1.nw minuit.nw postlude.nw

if DISTRIBUTION
PDFS = circe1.pdf 
else
PDFS = 
endif

PICTURES_PDF = \
    figures1.pdf fit11.pdf fit12.pdf \
    fit21.pdf fit22.pdf fit13.pdf \
    fit23.pdf fit15.pdf fit25.pdf \
    dist78.pdf

LATEX_STYLES = \
    noweb.sty thohacks.sty thopp.sty

TEX_FLAGS = "$$TEXINPUTS:$(top_srcdir)/circe1/share/doc"
MP_FLAGS  = "$$MPINPUTS:$(top_srcdir)/circe1/share/doc"

MP4_FILES = \
    circemacs.mp4 dist.mp4 fit.mp4 graphx.mp

EXTRA_DIST = \
  tex-comments.sh \
  $(PICTURES_PDF) \
  $(LATEX_STYLES) \
  $(MP4_FILES)

dist_doc_DATA = $(PDFS)

if NOWEB_AVAILABLE
pdf-local: circe1.pdf 
endif

VPATH = $(srcdir):$(top_builddir)/circe1/src:$(top_srcdir)/circe1/src

if NOWEB_AVAILABLE

circe1.tex: $(WEBS)
	@if $(AM_V_P); then :; else echo " NOWEAVE   " $@; fi
	$(AM_V_at)$(NOWEAVE) -filter ./tex-comments -delay -index \
	    `for i in  $^; do case $$i in *.nw) echo $$i;; esac done` \
	  > $@

circe1.tex: tex-comments

endif NOWEB_AVAILABLE

.mp4.mp: circemacs.mp4
	@if $(AM_V_P); then :; else echo " M4        " $@; fi
	$(AM_V_at)$(M4) -I$(top_srcdir)/circe1/share/doc $< >$@

########################################################################
# Old targets for the fit plots
########################################################################

fit.mp dist.mp: circemacs.mp4

fit.11: fit.mp
	TEX=$(LATEX) TEXINPUTS=$(TEX_FLAGS) $(MPOST) $<

########################################################################
# The distribution plots for the write-up
########################################################################

dist.1: dist.mp \
  de-sband.dat dg-sband.dat de-tesla.dat dg-tesla.dat \
  de-xband.dat dg-xband.dat de-sbandt.dat dg-sbandt.dat \
  de-teslat.dat dg-teslat.dat de-xbandt.dat dg-xbandt.dat \
  de-tesla3.dat dg-tesla3.dat de-tesla8.dat dg-tesla8.dat \
  de-sband-ee.dat dg-sband-ee.dat de-tesla-ee.dat dg-tesla-ee.dat \
  de-xband-ee.dat dg-xband-ee.dat de-sbandt-ee.dat dg-sbandt-ee.dat \
  de-teslat-ee.dat dg-teslat-ee.dat de-xbandt-ee.dat dg-xbandt-ee.dat
	@if $(AM_V_P); then TEX=$(LATEX) TEXINPUTS=$(TEX_FLAGS) \
	MPINPUTS=$(MP_FLAGS) $(MPOST) $<; else \
	echo " METAPOST  " $@; TEX=$(LATEX) TEXINPUTS=$(TEX_FLAGS) \
	MPINPUTS=$(MP_FLAGS) $(MPOST) $< >/dev/null; fi

########################################################################
# Generating CIRCE data
########################################################################

de-sband.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 11 11 500.0 1 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

dg-sband.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 22 11 500.0 1 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

de-tesla.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 11 11 500.0 2 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

dg-tesla.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 22 11 500.0 2 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

de-xband.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 11 11 500.0 3 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

dg-xband.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 22 11 500.0 3 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

de-sbandt.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 11 11 1000.0 1 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

dg-sbandt.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 22 11 1000.0 1 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

de-teslat.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 11 11 1000.0 2 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

dg-teslat.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 22 11 1000.0 2 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

de-xbandt.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 11 11 1000.0 3 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

dg-xbandt.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 22 11 1000.0 3 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

de-tesla3.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 11 11 350.0 2 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

dg-tesla3.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 22 11 350.0 2 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

de-tesla8.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 11 11 800.0 2 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

dg-tesla8.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 22 11 800.0 2 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

de-sband-ee.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 11 11 500.0 4 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

dg-sband-ee.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 22 11 500.0 4 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

de-tesla-ee.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 11 11 500.0 5 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

dg-tesla-ee.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 22 11 500.0 5 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

de-xband-ee.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 11 11 500.0 6 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

dg-xband-ee.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 22 11 500.0 6 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

de-sbandt-ee.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 11 11 1000.0 4 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

dg-sbandt-ee.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 22 11 1000.0 4 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

de-teslat-ee.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 11 11 1000.0 5 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

dg-teslat-ee.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 22 11 1000.0 5 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

de-xbandt-ee.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 11 11 1000.0 6 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

dg-xbandt-ee.dat: $(top_builddir)/circe1/tools/circe1_plot
	@if $(AM_V_P); then :; else echo " PLOT      " $@; fi
	$(AM_V_at)echo 1 0.001 1.0 50 -1.0 22 11 1000.0 6 $(CIRCE1_VER) $(CIRCE1_REV) | $(CPLOT)

CPLOT = $(top_builddir)/circe1/tools/circe1_plot | $(GREP) -v circe1: | $(SED) 's/E/e/g' > $@

tex-comments: tex-comments.sh
	cp $< $@
	chmod +x $@

SUFFIXES = .tex .pdf

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

if DISTRIBUTION
if SUPP_PDF_AVAILABLE
if PDFLATEX_AVAILABLE
if !NOWEB_AVAILABLE
circe1.pdf: $(PICTURES_PDF) dist.1
else NOWEB_AVAILABLE
circe1.pdf: $(PICTURES_PDF) dist.1 circe1.tex
endif NOWEB_AVAILABLE
endif PDFLATEX_AVAILABLE
endif SUPP_PDF_AVAILABLE
endif DISTRIBUTION


########################################################################

## Cleanup tasks
mostlyclean-latex:
	-rm -f *.data *.mpx *.[1-9] *.t[1-9] circe*.mp preview*.mp \
		circe1.tex *.out *.log *.aux *.idx *.ilg *.ind *.toc \
		$(PICTURES_PDF) tex-comments circe1_plot *.dat dist.mp \
		fit.mp dist.11 dist.12 dist.13 dist.14 
	-test "$(srcdir)" != "." && rm -f $(PDFS) 
clean-latex:
maintainer-clean-latex: 
	-rm -f $(PDFS)
if NOWEB_AVAILABLE
mostlyclean-circe1:
	-test "$(srcdir)" != "." && rm -f $(PDFS)
maintainer-clean-circe1:
else
mostlyclean-circe1:
maintainer-clean-circe1:
endif
.PHONY: mostlyclean-latex clean-latex maintainer-clean-latex
.PHONY: mostlyclean-circe1 maintainer-clean-circe1

## Remove backup files
maintainer-clean-backup:
	-rm -f *~
.PHONY: maintainer-clean-backup

## Register additional clean targets
mostlyclean-local: mostlyclean-latex mostlyclean-circe1
clean-local: clean-latex 
maintainer-clean-local: maintainer-clean-backup \
	maintainer-clean-circe1 maintainer-clean-latex

if !DISTRIBUTION
install-data-hook: 
	-$(INSTALL) -m 644 circe1.pdf $(DESTDIR)$(datarootdir)/doc/circe1

uninstall-hook: 
	-rm -f $(DESTDIR)/$(datarootdir)/doc/circe1/circe1.pdf 
endif

########################################################################
## The End.
########################################################################

