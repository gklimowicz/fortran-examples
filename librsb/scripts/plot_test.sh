#!/bin/bash
#
# Copyright (C) 2020-2020 Michele Martone
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

# This script is intended the librsb developer usage.
set -e
set -x
( ! test -d tmpdir ) || { echo "Please remove or rename tmpdir! "; exit 1; } 

( ! gnuplot -c plot.gpi --non-existing-option ) # fail on error
if test -n "$DISPLAY" ; then
	gnuplot -c plot.gpi --no-plot # this relies on X11
fi
gnuplot -c plot.gpi --dumb
yes | gnuplot -c plot.gpi --pause
gnuplot -c plot.gpi 
rm -f non-existent-file
test `gnuplot -c plot.gpi    --png -H --datafile non-existent-file 2>&1  | wc -l ` = 1
gnuplot -c plot.gpi --logfile plot_sample.log --no-plot --all --mtxgrep --idx 0  2>&1 | grep Found.matrices

( ! gnuplot -c plot.gpi --datafile plot_sample.dat --no-plot --all --idx 2 ) # idx=2 plot_sample.dat's two blocks (0 1)

TMPDIR=`mktemp -d  /dev/shm/gnuplot.tmpdir.XXXXX`
test -d "$TMPDIR" && trap "rm -fR $TMPDIR; rm tmpdir; true; " EXIT
ln -s $TMPDIR tmpdir

rm -f tmpdir/temp*png
cp plot_sample.log temp.log
gnuplot -c plot.gpi --png --logfile temp.log  --outdir tmpdir --all --idx 0 --boxes
test -f tmpdir/temp_RSBBEST-MFLOPS_for__all.png
test -f tmpdir/temp_RSBBEST-MFLOPS_for__all_boxes.png
test -f tmpdir/temp_AT-OPTIME_for__all.png
test -f tmpdir/temp_AT-OPTIME_for__all_boxes.png
rm   -f tmpdir/temp_*png
rm -f temp.log

rm -f tmpdir/temp*png
cp plot_sample.log temp.log
gnuplot -c plot.gpi --png --logfile temp.log  --outdir tmpdir --all --idx '0:1' # or '0 1' (still an interval)
test -f tmpdir/temp_0_RSBBEST-MFLOPS_for__all.png
! test -f tmpdir/temp_0_RSBBEST-MFLOPS_for__all_boxes.png
test -f tmpdir/temp_0_AT-OPTIME_for__all.png
! test -f tmpdir/temp_0_AT-OPTIME_for__all_boxes.png
test -f tmpdir/temp_1_RSBBEST-MFLOPS_for__all.png
! test -f tmpdir/temp_1_RSBBEST-MFLOPS_for__all_boxes.png
test -f tmpdir/temp_1_AT-OPTIME_for__all.png
! test -f tmpdir/temp_1_AT-OPTIME_for__all_boxes.png
rm   -f tmpdir/temp_*png
rm -f temp.log

rm -f tmpdir/temp*png
cp plot_sample.log temp.log
gnuplot -c plot.gpi --png --logfile temp.log  --outdir tmpdir --all --idx ':' # or '*'
test -f tmpdir/temp_0_RSBBEST-MFLOPS_for__all.png
! test -f tmpdir/temp_0_RSBBEST-MFLOPS_for__all_boxes.png
test -f tmpdir/temp_0_AT-OPTIME_for__all.png
! test -f tmpdir/temp_0_AT-OPTIME_for__all_boxes.png
test -f tmpdir/temp_1_RSBBEST-MFLOPS_for__all.png
! test -f tmpdir/temp_1_RSBBEST-MFLOPS_for__all_boxes.png
test -f tmpdir/temp_1_AT-OPTIME_for__all.png
! test -f tmpdir/temp_1_AT-OPTIME_for__all_boxes.png
rm   -f tmpdir/temp_*png
rm -f temp.log

rm -f tmpdir/temp*png
cp plot_sample.dat temp.dat
gnuplot -c plot.gpi --png --datafile temp.dat --outdir tmpdir --idx 0 --all
test -f tmpdir/temp_RSBBEST-MFLOPS_for__all.png
! test -f tmpdir/temp_RSBBEST-MFLOPS_for__all_boxes.png
test -f tmpdir/temp_AT-OPTIME_for__all.png
! test -f tmpdir/temp_AT-OPTIME_for__all_boxes.png
rm -f temp.dat

rm -f tmpdir/temp*.png
cp plot_sample.dat temp.dat
gnuplot -c plot.gpi --png --datafile temp.dat --mtx dense-10000x10000-100000000nz --outdir tmpdir --idx 0
test ! -f tmpdir/temp_AT-OPTIME_for_MTX_COPY1_10000x10000-100000000nz.png
test ! -f tmpdir/temp_AT-OPTIME_for_MTX_COPY2_10000x10000-100000000nz.png
test   -f tmpdir/temp_AT-OPTIME_for_MTX_dense-10000x10000-100000000nz.png
test ! -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_COPY1_10000x10000-100000000nz.png
test ! -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_COPY2_10000x10000-100000000nz.png
test   -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_dense-10000x10000-100000000nz.png
rm -f temp.dat

rm -f tmpdir/temp*.fig
cp plot_sample.dat temp.dat
gnuplot -c plot.gpi --fig --datafile temp.dat --mtx dense-10000x10000-100000000nz --outdir tmpdir --idx 0
test ! -f tmpdir/temp_AT-OPTIME_for_MTX_COPY1_10000x10000-100000000nz.fig
test ! -f tmpdir/temp_AT-OPTIME_for_MTX_COPY2_10000x10000-100000000nz.fig
test   -f tmpdir/temp_AT-OPTIME_for_MTX_dense-10000x10000-100000000nz.fig
test ! -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_COPY1_10000x10000-100000000nz.fig
test ! -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_COPY2_10000x10000-100000000nz.fig
test   -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_dense-10000x10000-100000000nz.fig
rm -f temp.dat

rm -f tmpdir/temp*.eps
cp plot_sample.dat temp.dat
gnuplot -c plot.gpi --eps --datafile temp.dat --mtx dense-10000x10000-100000000nz --outdir tmpdir --idx 0
test ! -f tmpdir/temp_AT-OPTIME_for_MTX_COPY1_10000x10000-100000000nz.eps
test ! -f tmpdir/temp_AT-OPTIME_for_MTX_COPY2_10000x10000-100000000nz.eps
test   -f tmpdir/temp_AT-OPTIME_for_MTX_dense-10000x10000-100000000nz.eps
test ! -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_COPY1_10000x10000-100000000nz.eps
test ! -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_COPY2_10000x10000-100000000nz.eps
test   -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_dense-10000x10000-100000000nz.eps
rm -f temp.dat

rm -f tmpdir/temp.*jpg
cp plot_sample.dat temp.dat
gnuplot -c plot.gpi --jpg --datafile temp.dat --mtx dense-10000x10000-100000000nz --outdir tmpdir --idx 0
test ! -f tmpdir/temp_AT-OPTIME_for_MTX_COPY1_10000x10000-100000000nz.jpg
test ! -f tmpdir/temp_AT-OPTIME_for_MTX_COPY2_10000x10000-100000000nz.jpg
test   -f tmpdir/temp_AT-OPTIME_for_MTX_dense-10000x10000-100000000nz.jpg
test ! -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_COPY1_10000x10000-100000000nz.jpg
test ! -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_COPY2_10000x10000-100000000nz.jpg
test   -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_dense-10000x10000-100000000nz.jpg
rm -f temp.dat

rm -f tmpdir/temp.*pdf
cp plot_sample.dat temp.dat
gnuplot -c plot.gpi --pdf --datafile temp.dat --mtx dense-10000x10000-100000000nz --outdir tmpdir --idx 0
test ! -f tmpdir/temp_AT-OPTIME_for_MTX_COPY1_10000x10000-100000000nz.pdf
test ! -f tmpdir/temp_AT-OPTIME_for_MTX_COPY2_10000x10000-100000000nz.pdf
test   -f tmpdir/temp_AT-OPTIME_for_MTX_dense-10000x10000-100000000nz.pdf
test ! -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_COPY1_10000x10000-100000000nz.pdf
test ! -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_COPY2_10000x10000-100000000nz.pdf
test   -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_dense-10000x10000-100000000nz.pdf
rm -f temp.dat

rm -f tmpdir/temp*.svg
cp plot_sample.dat temp.dat
gnuplot -c plot.gpi --svg --datafile temp.dat --mtx dense-10000x10000-100000000nz --outdir tmpdir --idx 0
test ! -f tmpdir/temp_AT-OPTIME_for_MTX_COPY1_10000x10000-100000000nz.svg
test ! -f tmpdir/temp_AT-OPTIME_for_MTX_COPY2_10000x10000-100000000nz.svg
test   -f tmpdir/temp_AT-OPTIME_for_MTX_dense-10000x10000-100000000nz.svg
test ! -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_COPY1_10000x10000-100000000nz.svg
test ! -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_COPY2_10000x10000-100000000nz.svg
test   -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_dense-10000x10000-100000000nz.svg
rm -f temp.dat

rm -f tmpdir/temp*.tex
cp plot_sample.dat temp.dat
gnuplot -c plot.gpi --tex            temp.dat --mtx dense-10000x10000-100000000nz --outdir tmpdir --idx 0
test ! -f tmpdir/temp_AT-OPTIME_for_MTX_COPY1_10000x10000-100000000nz.tex
test ! -f tmpdir/temp_AT-OPTIME_for_MTX_COPY2_10000x10000-100000000nz.tex
test   -f tmpdir/temp_AT-OPTIME_for_MTX_dense-10000x10000-100000000nz.tex
test ! -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_COPY1_10000x10000-100000000nz.tex
test ! -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_COPY2_10000x10000-100000000nz.tex
test   -f tmpdir/temp_RSBBEST-MFLOPS_for_MTX_dense-10000x10000-100000000nz.tex
rm -f temp.dat

#cp plot_sample.dat plot_sample.dat
rm   -f plot_sample_*.png
gnuplot -c plot.gpi --png --datfile plot_sample.dat --idx 0
test -f plot_sample_AT-OPTIME_for_MTX_COPY1_10000x10000-100000000nz.png
test -f plot_sample_AT-OPTIME_for_MTX_COPY2_10000x10000-100000000nz.png
test -f plot_sample_AT-OPTIME_for_MTX_dense-10000x10000-100000000nz.png
test -f plot_sample_AT-OPTIME_for_NT_24.png
test -f plot_sample_AT-OPTIME_for_NT_32.png
test -f plot_sample_AT-OPTIME_for_TRANS_N.png
test -f plot_sample_AT-OPTIME_for_TRANS_T.png
test -f plot_sample_AT-OPTIME_for_TYPE_C.png
test -f plot_sample_AT-OPTIME_for_TYPE_D.png
test -f plot_sample_AT-OPTIME_for_TYPE_S.png
test -f plot_sample_AT-OPTIME_for_TYPE_Z.png
test -f plot_sample_RSBBEST-MFLOPS_for_MTX_COPY1_10000x10000-100000000nz.png
test -f plot_sample_RSBBEST-MFLOPS_for_MTX_COPY2_10000x10000-100000000nz.png
test -f plot_sample_RSBBEST-MFLOPS_for_MTX_dense-10000x10000-100000000nz.png
test -f plot_sample_RSBBEST-MFLOPS_for_NT_24.png
test -f plot_sample_RSBBEST-MFLOPS_for_NT_32.png
test -f plot_sample_RSBBEST-MFLOPS_for_TRANS_N.png
test -f plot_sample_RSBBEST-MFLOPS_for_TRANS_T.png
test -f plot_sample_RSBBEST-MFLOPS_for_TYPE_C.png
test -f plot_sample_RSBBEST-MFLOPS_for_TYPE_D.png
test -f plot_sample_RSBBEST-MFLOPS_for_TYPE_S.png
test -f plot_sample_RSBBEST-MFLOPS_for_TYPE_Z.png
rm -f plot_sample*png
echo 'All tests passed.'
