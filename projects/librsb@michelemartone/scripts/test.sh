#!/bin/bash
#
# Copyright (C) 2008-2021 Michele Martone
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

# Script for library consistence checking.

set -e
set -x
set -o pipefail
shopt -s expand_aliases

if test x"${srcdir}" = x ; then srcdir=. ; fi
ULIMIT_S=10000
echo "Invoking ulimit -s ${ULIMIT_S}"
ulimit -s ${ULIMIT_S} 


value="L2:4/64/512K,L1:8/64/24K";
test `RSB_USER_SET_MEM_HIERARCHY_INFO=$value ./rsbench -C | grep RSB_USER_SET_MEM_HIERARCHY_INFO` == \
      RSB_USER_SET_MEM_HIERARCHY_INFO:$value
unset value;

test `RSB_USER_SET_MEM_HIERARCHY_INFO=value  ./rsbench -C | grep RSB_USER_SET_MEM_HIERARCHY_INFO` == \
      RSB_USER_SET_MEM_HIERARCHY_INFO:value

# check rsb__get_lnc_size() output:
O=`( set +o pipefail; OMP_NUM_THREADS=1 RSB_USER_SET_MEM_HIERARCHY_INFO="L2:4/64/512K,L1:8/64/24K" ./rsbench -I | less | head -n 10 ) | grep '^L. size:' `
echo "${O}" 
test "${O}" == "L1 size: 24576 
L2 size: 524288 "

QUICK='--no-flush-cache-in-iterations --no-flush-cache-around-loop --want-no-memory-benchmark'

gen_liminal()
{
	ms=`$B --configuration | grep RSB_MAX_MATRIX_DIM |sed "s/^RSB_MAX_MATRIX_DIM://g"` ; 
	[ -z "$ms" ] && exit 1
	r=$ms
	c=$ms
	n=4
cat << EOF > $M
%%MatrixMarket matrix coordinate real general
% this matrix is sized as the limits of your architecture allow
$r $c $n
1 1
1 $c
$r 1
$r $c
EOF
}
nulldev=/dev/null
alias 2>&1 > $nulldev || exit 1

alias fail='exit 1'
#alias fail='( echo "[!]" ; exit 1 ; )' # seems to not fail..
e='echo'
x='echo "--"; echo'
strict=1
d=`pwd`
B=$d/rsbench
M=$d/test.mtx
n=100
pt=0
ft=0
st=0
fc=""
$x || exit -1

make rsbench || fail	# we want our test programs

test `$B -C  | grep row.unrolls |sed "s/^row unrolls://g;s/ /,/g"` = 1
test `$B -C  | grep column.unrolls |sed "s/^column unrolls://g;s/ /,/g"` = 1
test `$B -C  | grep ^types.count: | sed s/types.count://g` -gt 0
test `$B -C  | grep ^types.count: | sed s/types.count://g` -lt 8

$B -oa -Ob ${QUICK} --want-no-autotune --notranspose -f ${srcdir}/us.mtx --expand-symmetry --verbose  | grep us.mtx.*36
$B -oa -Ob ${QUICK} --matrix-sample-pcnt   1 --want-no-autotune --no-transpose -f ${srcdir}/A.mtx | grep sing.only.1. # alias to --notranspose
$B -oa -Ob ${QUICK} --matrix-sample-pcnt   1 --want-no-autotune --notranspose -f ${srcdir}/A.mtx | grep sing.only.1.
$B -oa -Ob ${QUICK} --matrix-sample-pcnt  50 --want-no-autotune --notranspose -f ${srcdir}/A.mtx | grep sing.only.3.
$B -oa -Ob ${QUICK} --matrix-sample-pcnt  99 --want-no-autotune --notranspose -f ${srcdir}/A.mtx | grep sing.only.5.

# FIXME : compact all of these flags in some way..
if true ; then 

$x "$B --help" # -h
    $B --help || fail # 

$x "$B -oa -Ob -h"
    $B -oa -Ob -h || fail # 

#$x "$B -M"
#    $B -M || fail # 

$x "$B --hardware-counters" # -H
    $B --hardware-counters || fail # 

$x "$B --configuration" # -C
    $B --configuration || fail # 

pdm=${srcdir}/pd.mtx

$x "$B --guess-blocking $pdm" # -G
    $B --guess-blocking $pdm || fail # 

$x "$B -oa -Ob  ${QUICK} -f $pdm --matrix-dump-graph $pdm.dot --matrix-dump-internals"
$B -oa -Ob  ${QUICK} -f $pdm --matrix-dump-graph $pdm.dot --matrix-dump-internals || fail # 

$x "$B -ot -Ob --lower 3"
$B -ot -Ob --lower 3 || fail # 

$x "$B --matrix-print $pdm" # -P
    $B --matrix-print $pdm || fail # 

$x "$B -I"
    $B -I || fail # system information dumpout

bmfn=test.mtx.rsb
for deff in "-R -Fbo -qH" "-R -Fbo" ; do
for detr in "--lower" "--dense" ; do
if $B --configuration | grep 'XDR.*off' ; then
	st=$((st+1))
else
	$x "$B -oa -Ob ${QUICK} $detr=10 -w $bmfn $deff"
	    $B -oa -Ob ${QUICK} $detr=10 -w $bmfn  || fail # binary I/O test
	
	$x "$B -oa -Ob ${QUICK} -b $bmfn $deff"
	    $B -oa -Ob ${QUICK} -b $bmfn  || fail # binary I/O test
fi
done
done

$x "$B -OR"
    $B -OR || fail # a no op

$x "$B -Ot -b"
    $B -Ot -b || fail # a no op

$x "$B --bench -e -f test.mtx"
    $B --bench -e -f test.mtx || fail # obsolete; here for coverage

$x "$B  -oa -Ob ${QUICK} --gen-diag 10000  --beta 2 --implicit-diagonal"
    $B  -oa -Ob ${QUICK} --gen-diag 10000  --beta 2 --implicit-diagonal || fail

if test x"$RSB_SHORT_TEST_SH" = x1; then
	export RSB_BENCHMARK_MIN_SECONDS=0.01
else
	export RSB_BENCHMARK_MIN_SECONDS=0.1
fi
$x "$B -Or" # reference
    $B -Or || fail

$x "$B --matrix-ls-latex ${srcdir}/A.mtx"
    $B --matrix-ls-latex ${srcdir}/A.mtx || fail #
diff <( $B --matrix-ls-latex ${srcdir}/A.mtx ; ) - <<- EOF
\\begin{table}[]\\begin{footnotesize}\\begin{center} \\begin{tabular}{lllll}\\hline
matrix & rows & columns & nnz & nnz/row \\\\\\hline
A.mtx & 3 & 3 & 6 & 2\\\\%%symm
\\hline \\end{tabular} \\caption{Caption.}\\label{testbed_matrices}\\end{center}\\end{footnotesize}\\end{table}
EOF

$x "$B -Oc -f $M"
RSB_BENCHMARK_MIN_SECONDS=0.0001 $B -Oc -f $M || fail # fuller testing

	$B --plot-matrix -aRzd -f $pdm || fail
fi

# The following are here for (rather flimsy) testing/coverage purposes:
$x "$B -oa -Ob ${QUICK} -f ${srcdir}/A.mtx -Fo  --z-sorted-coo"
$B -oa -Ob ${QUICK} -f ${srcdir}/A.mtx -Fo  --z-sorted-coo || fail
$x "$B -oa -Ob ${QUICK} --lower 4 --want-no-recursive --ilu0"
$B -oa -Ob ${QUICK} --lower 4 --want-no-recursive --ilu0 || fail
$x "$B -oa -Ob ${QUICK} --lower 4 -K"
$B -oa -Ob ${QUICK} --lower 4 -K || fail
$x "$B -oa -Ob ${QUICK} --dense 10 --nrhs 1,2 --incy 1,2 --nrhs 1,2 -K"
$B -oa -Ob ${QUICK} --dense 10 --nrhs 1,2 --incy 1,2 --nrhs 1,2 -K

if make sbtf -a -f sbtf ; then ./sbtf || fail ; fi # make sbtf gives no error code if sbtf conditionally out.
if make sbtc -a -f sbtc ; then ./sbtc || fail ; fi
if make ot   -a -f   ot ; then ./ot   || fail ; fi

for o in "-b 9" "-n 10%" "-n 2%"; # bandwidth and overall density
do
for n in 1 2 3 4  10 100 200;
do
case $n in 
	-1)
	gen_liminal ;; # in one case we test for a limit sized matrix
	1|2|3|4)
	$x "$B --generate-matrix -r $n -c $n -b $n  > $M" # -g
	   $B --generate-matrix -r $n -c $n -b $((n-1))  > $M || fail
	;;
	*)
	$x "$B --generate-matrix -r $n -c $n $o  > $M"
	   $B --generate-matrix -r $n -c $n $o  > $M || fail
esac

# provide a small matrix when matrix creation was disabled at configure time
$B --configuration | grep RSB_IOLEVEL:0 && cp $pdm $M 

sep="              *"

for f in `$B --configuration | grep format.switches |sed "s/^format switches://g"` ; # for every supported format
do

# various formats testing
for a in "" "-R";	# with RSB format
do
	# matrix dumpout
	$x "$B -Od -f $M -T : -F $f"
	   $B -Od -f $M -T : -F $f > $nulldev || fail # needs -f matrixname

#	$x "$B -oa -Ob ${QUICK} -f $M -d -F $f -t 1 $a"
#	   $B -oa -Ob ${QUICK} -f $M -d -F $f -t 1 $a || fail
	c="$?"
#	$x "$B -oa -Ob ${QUICK} -f $M -d -F $f -t 1 $a"
#	es="$B -oa -Ob ${QUICK} -f $M -d -F $f -t 1 $a"

	$x "$B -oa -Ob ${QUICK} -f $M  -F $f -t 1 $a"
	   $B -oa -Ob ${QUICK} -f $M  -F $f -t 1 $a || fail
	c="$?"
	$x "$B -oa -Ob ${QUICK} -f $M  -F $f -t 1 $a"
	es="$B -oa -Ob ${QUICK} -f $M  -F $f -t 1 $a"

	if test x"$c" = x"0" 
	then
		pt=$((pt+1))
	else
		ft=$((ft+1))
		fc="$fc\n$es"
		if test x"$strict" = x"1" ; then
			# FIXME
			exit -1
		fi
	fi
done
done
done
done

$x "$B -oa -Ob ${QUICK} -R --dense 100 --write-performance-record=test.rpr"
    $B -oa -Ob ${QUICK} -R --dense 100 --write-performance-record=test.rpr
$x "$B                      --read-performance-record=test.rpr"
    $B                      --read-performance-record=test.rpr

$x "$B -oa -Ob ${QUICK} -R --write-performance-record=test.rpr ${pdm} non-existing.mtx"
    $B -oa -Ob ${QUICK} -R --write-performance-record=test.rpr ${pdm} non-existing.mtx

for RSB_PR_MULTIDUMP in 0 1 2 3; do
	RSB_PR_MULTIDUMP=$RSB_PR_MULTIDUMP $B --read-performance-record  \
		test.rpr  test.rpr
done

rm test.rpr

# test for the RSB_FAF_CHKGSS feature
rm -f non-existing.mtx
rm -f non-existing.mtx.gz
cp -pv ${pdm} non-existing.mtx.gz
if grep RSB_WANT_ZLIB_SUPPORT.1 rsb-config.h ; then
	$x "$B -oa -Ob -R non-existing.mtx"
	    { $B -oa -Ob -R non-existing.mtx | grep assuming.you; } || fail;
fi

if ! test x"$RSB_SHORT_TEST_SH" = x1; then
$x "$B --blas-testing" # -B
    $B --blas-testing || fail # 
fi

$e "passed  tests : $pt"
$e "failed  tests : $ft"
$e "skipped tests : $st"

if test "$ft" != "0" ; then fail ; fi

echo $fc
