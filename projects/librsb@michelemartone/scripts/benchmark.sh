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
set -e
set -x
which rsbench
TRF=`mktemp librsb.reference.XXX.tmp`
touch ${TRF}
if test $# = 0; then
	./scripts/matrices_get.sh
fi
rsbench -oa -Ob --bench `if test $# = 0; then echo *.mtx ; else echo ${@}; fi` \
	--verbose --verbose
BRF=`find -name '*.rpr' -newer ${TRF}`
rm ${TRF}
rsbench --read-performance-record ${BRF} > ${BRF/rpr/txt}
ls -ltr ${BRF/rpr/txt}
RSB_PR_WLTC=2 RSB_PR_SR=0 \
 rsbench --read-performance-record ${BRF} > ${BRF/rpr/tex}
which pdflatex || exit 0
which kpsepath || exit 0 # to check LaTeX packages
find `kpsepath tex | sed 's/!!//g;s/:/\n/g;'` -name sciposter.cls || exit 0 # need sciposter class, usually in texlive-science
find `kpsepath tex | sed 's/!!//g;s/:/\n/g;'` -name translator.sty || exit 0 # need sciposter class, usually in texlive-latex-recommended
pdflatex -interaction=batchmode -file-line-error ${BRF/rpr/tex}
which gnuplot
RSB_PRD_STYLE_PLT_PFN=${BRF/rpr/} RSB_PRD_STYLE_PLT_FMT=1 RSB_PR_SR=2 \
 rsbench --read-performance-record ${BRF} > ${BRF/rpr/gnu}
ls -ltr ${BRF/rpr/gnu}
gnuplot ${BRF/rpr/gnu}
RSB_PRD_STYLE_PLT_PFN=${BRF/rpr/} RSB_PRD_STYLE_PLT_FMT=  RSB_PR_SR=2 \
 rsbench --read-performance-record ${BRF} > ${BRF/rpr/gnu}
gnuplot ${BRF/rpr/gnu}
rm ${BRF/rpr/}{aux,log,out,gnu}
ls -ltr ${BRF/.rpr/}*.png
ls -ltr ${BRF/.rpr/}*.eps
rm ${BRF/.rpr/}*{.png,.eps}
ls -ltr ${BRF} ${BRF/rpr/tex} ${BRF/rpr/pdf} ${BRF/rpr/txt}
rm      ${BRF} ${BRF/rpr/tex} ${BRF/rpr/pdf} ${BRF/rpr/txt}
exit
