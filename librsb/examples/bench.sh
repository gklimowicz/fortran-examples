#!/bin/bash
#
# Copyright (C) 2008-2022 Michele Martone
# 
# This file is part of librsb.
# 
# librsb is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# 
# librsb is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with librsb; see the file COPYING.
# If not, see <http://www.gnu.org/licenses/>.

# Benchmark rsbench and postprocess results.
# \ingroup rsb_doc_examples
# @file
# @author Michele Martone
# @brief Benchmark invocation from shell script.
#
set -e
set -x
which rsbench
BRF=test.rpr
# invoke rsbench and produce a performance record, using all types and one thread
rsbench -oa -Ob --bench --lower 100 --as-symmetric \
	--types ':' -n 1 \
	--notranspose --compare-competitors \
	--verbose --verbose \
	--write-performance-record=${BRF}

# examine tuning renderings (produced by --verbose --verbose)
ls -ltr ${BRF/.rpr/}-tuning*

# convert the performance record to text form
rsbench --read-performance-record ${BRF} > ${BRF/rpr/txt}
ls -ltr ${BRF/rpr/txt}

# convert the performance record to LaTeX table document form
RSB_PR_WLTC=2 RSB_PR_SR=0 \
 rsbench --read-performance-record ${BRF} > ${BRF/rpr/tex}
which latex || exit 0 # to compile LaTeX document
which kpsepath || exit 0 # to check LaTeX packages
find `kpsepath tex | sed 's/!!//g;s/:/\n/g;'` -name sciposter.cls || exit 0 # need sciposter class, usually in texlive-science
find `kpsepath tex | sed 's/!!//g;s/:/\n/g;'` -name translator.sty || exit 0 # need sciposter class, usually in texlive-latex-recommended

# convert the LaTeX table into a DVI (may as well use pdflatex for PDF)
latex -interaction=batchmode -file-line-error ${BRF/rpr/tex}

# convert the performance record to GNUPLOT plots
which gnuplot || exit 0
RSB_PRD_STYLE_PLT_PFN=${BRF/rpr/} RSB_PRD_STYLE_PLT_FMT=1 RSB_PR_SR=2 \
 rsbench --read-performance-record ${BRF} > ${BRF/rpr/gnu}

# convert the GNUPLOT plots into PDF
ls -ltr ${BRF/rpr/gnu}
gnuplot ${BRF/rpr/gnu}

# convert the performance record to GNUPLOT plots, different way
RSB_PRD_STYLE_PLT_PFN=${BRF/rpr/} RSB_PRD_STYLE_PLT_FMT=  RSB_PR_SR=2 \
 rsbench --read-performance-record ${BRF} > ${BRF/rpr/gnu}
gnuplot ${BRF/rpr/gnu}
ls -ltr ${BRF/.rpr/}*.png
ls -ltr ${BRF/.rpr/}*.eps
ls -ltr ${BRF} ${BRF/rpr/tex} ${BRF/rpr/dvi} ${BRF/rpr/txt}

# clean up
rm ${BRF/rpr/}{aux,log,out,gnu}
rm ${BRF/.rpr/}*{.png,.eps}
rm      ${BRF} ${BRF/rpr/tex} ${BRF/rpr/dvi} ${BRF/rpr/txt}
exit

