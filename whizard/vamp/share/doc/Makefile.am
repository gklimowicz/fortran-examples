# Makefile.am --
########################################################################
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


WEBS = \
  prelude.nw divisions.nw vamp.nw vampi.nw \
  vamp_test.nw vamp_test0.nw application.nw \
  vamp_kinds.nw constants.nw exceptions.nw \
  tao_random_numbers.nw specfun.nw vamp_stat.nw histograms.nw \
  utils.nw linalg.nw products.nw kinematics.nw coordinates.nw \
  mpi90.nw postlude.nw

if DISTRIBUTION
PDFS = vamp.pdf preview.pdf preview2.pdf
else
PDFS = 
endif

LATEX_STYLES = \
    feynmp.sty feynmp.mp \
    noweb.sty emp.sty flex.cls thohacks.sty thophys.sty

TEX_FLAGS = "$$TEXINPUTS:$(top_srcdir)/vamp/share/doc"
MP_FLAGS  = "$$MPINPUTS:$(top_srcdir)/vamp/share/doc"

EXTRA_DIST = \
  tex-comments.sh \
  vegas.d vamp.d \
  $(LATEX_STYLES)

dist_doc_DATA = $(PDFS)

if NOWEB_AVAILABLE
pdf-local: vamp.pdf preview.pdf preview2.pdf
else
pdf-local: preview.pdf preview2.pdf
endif

VPATH = $(srcdir):$(top_builddir)/vamp/src:$(top_srcdir)/vamp/src

if NOWEB_AVAILABLE

vamp.tex: $(WEBS)
	@if $(AM_V_P); then :; else echo " NOWEAVE   " $@; fi
	$(AM_V_at)$(NOWEAVE) -filter ./tex-comments -delay -index \
	    `for i in  $^; do case $$i in *.nw) echo $$i;; esac done` \
	  | $(CPIF) $@

vamp.tex: tex-comments

endif NOWEB_AVAILABLE

tex-comments: tex-comments.sh
	cp $< $@
	chmod +x $@

preview.pdf: vegas.data vamp.data

vegas.data: vegas.d
	cp $< $@

vamp.data: vamp.d
	cp $< $@

SUFFIXES = .tex .pdf

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
	@if $(AM_V_P); then test -r $*pics.mp && MPINPUTS=$(MP_FLAGS) $(MPOST) $*pics; else \
	   echo " METAPOST  " $*pics.mp; \
	   test -r $*pics.mp && MPINPUTS=$(MP_FLAGS) $(MPOST) $*pics > /dev/null; fi
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
	-rm -f *.data *.mpx *.[1-9] *.t[1-9] vamp*.mp preview*.mp \
	  *.out *.log *.aux *.idx *.ilg *.ind *.rcs *.toc \
	  tex-comments vamp.tex
	-test "$(srcdir)" != "." && rm -f vamp.pdf \
		preview.pdf preview2.pdf
clean-latex:
maintainer-clean-latex:
	-rm -f vamp.pdf preview.pdf preview2.pdf
if NOWEB_AVAILABLE
mostlyclean-vamp:
	-test "$(srcdir)" != "." && rm -f vamp.pdf \
		preview.pdf preview2.pdf
maintainer-clean-vamp:
else
mostlyclean-vamp:
maintainer-clean-vamp:
endif
.PHONY: mostlyclean-latex clean-latex maintainer-clean-latex
.PHONY: mostlyclean-vamp maintainer-clean-vamp

## Remove backup files
maintainer-clean-backup:
	-rm -f *~
.PHONY: maintainer-clean-backup

## Register additional clean targets
mostlyclean-local: mostlyclean-latex mostlyclean-vamp
clean-local: clean-latex
maintainer-clean-local: maintainer-clean-latex maintainer-clean-vamp \
	maintainer-clean-backup

if !DISTRIBUTION
install-data-hook: 
	-$(INSTALL) -m 644 vamp.pdf $(DESTDIR)$(datarootdir)/doc/vamp 
	-$(INSTALL) -m 644 preview.pdf $(DESTDIR)$(datarootdir)/doc/vamp 
	-$(INSTALL) -m 644 preview2.pdf $(DESTDIR)$(datarootdir)/doc/vamp 

uninstall-hook: 
	-rm -f $(DESTDIR)/$(datarootdir)/doc/vamp/vamp.pdf 
	-rm -f $(DESTDIR)/$(datarootdir)/doc/vamp/preview.pdf 
	-rm -f $(DESTDIR)/$(datarootdir)/doc/vamp/preview2.pdf
endif

########################################################################
## The End.
########################################################################

