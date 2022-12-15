#!/bin/bash
#
# Copyright (C) 2008-2015 Michele Martone
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

#
# This is a librsb developer tool, so don't expect documentation for it.
#
echo "#"
export RSBDB_COMPARE_SWITCH=${RSBDB_COMPARE_SWITCH:---compare-competitors}
export RSBDB_TIMES=${RSBDB_TIMES:-30}
export RSBDB_RSBENCH=${RSBDB_RSBENCH:-./rsbench}
export RSBDB_LABEL=${RSBDB_LABEL:-}
export RSBDB_LISTFILE=${RSBDB_LISTFILE:-}
#version=`grep '\<VERSION\>' config.h | sed 's/.*VERSION *//g;s/"//g;s/0\.//g'`
#pversion=`grep '\<PACKAGE_VERSION\>' config.h | sed 's/.*VERSION *//g;s/"//g;s/0\.//g'`
pversion=`./rsbench -C | grep version| sed 's/^.*version *: *//g;s/  */ /g' `
export RSBDB_VERSION=${RSBDB_VERSION:-$pversion}
#export RSBDB_MINDIMLOG=0
#export RSBDB_MAXDIMLOG=11
export RSBDB_MINDIMLOG=${RSBDB_MINDIMLOG:-5}
if cat /proc/cpuinfo | grep Xeon ; then
export RSBDB_MAXDIMLOG=${RSBDB_MAXDIMLOG:-12}
else
export RSBDB_MAXDIMLOG=${RSBDB_MAXDIMLOG:-11}
fi
echo "#"
env | grep ^RSBDB_
echo "#"
export RSBDB_FILESLIST=""
#
if test "$#" != 1 -o x"$1" != x" " ; then echo -e "Invocation example: \nRSBDB_RSBENCH=./rsbench $0 ' ' "; exit; fi
#
function proj(){  sed "s/\s\s*/ /g;s/^..//g;s/:/ /g"  | cut -f 1,2,5,9  -d  \  | sed "s/ / /" | grep PERFORMANCE| grep -v SERIAL | sed 's/\<PERFORMANCE/RSB_PERFORMANCE/g' ; }
function rproj() { cut -d \  -f 2,4 ; }
function plotstuff()
{
export RSBDB_FILESLIST="${RSBDB_FILESLIST} $1"
shopt -s expand_aliases
pentries=`cat "$1" | proj | cut -f 1 -d \   | sort | uniq`
nthreads=`cat "$1" | proj | cut -f 3 -d \   | sort | uniq`
pbins=`cat "$1" | proj `
#echo "$pbins"
#echo $pentries $nthreads
pf=${1/log/gp}
{
of=${1/log/eps}
echo "set term postscript color eps enhanced "
echo "set title \""${1}\\n${cinfo}"\" "
echo "set logscale x"
echo "set key left Left"
echo "set xtics autofreq 2"
echo "set xlabel \"matrix dimension\""
echo "set ylabel \"MFLOPS\""
echo "set output \"${of}\""
echo -en plot
for nt in $nthreads ; do
for ne in $pentries ; do
	ptitle=`echo $ne-$nt | sed s/_/-/g`
	echo -en ' "-" u 1:2 title "'$ptitle'" with linespoints ,'
done
done | sed 's/,$//g'
echo
for nt in $nthreads ; do
for ne in $pentries ; do
	FILTER="^$ne [[:graph:]][[:graph:]]* $nt "
	cat "$1" |  proj | grep "${FILTER}"  | rproj
	echo "#grep ${FILTER}:"
	echo e
	echo 
done
done
} > $pf
gnuplot $pf
}

# have:
# * dense SPMV
# * dense SPMV, with varying cache blocking
# * dense-spaced SPMV
# * dense SPMV, COO
# * dense SPMV, CSR
# * dense SPMV with incx
# * dense SPSV

# need:
# * dense SPSV, with varying cache blocking
# * types

CFLAGS=`./rsbench -I | grep CFLAGS| sed 's/CFLAGS *: *//g;s/  */ /g'  `
CC=`./rsbench -I | grep CC| sed 's/CC *: *//g;s/  */ /g'`
MKLINFO=`./rsbench   -oa -Ob   -R --lower 10 -qH --verbose  ${RSBDB_COMPARE_SWITCH} | grep '#%:MKL'`
export cinfo="`echo ${CC} ${CFLAGS} // ${MKLINFO} | sed 's/\(.\{80\}\)/\1\\\n/g'`"

if false ; then
nproc="`./rsbench -I | grep processors.online | sed 's/.*://' | sed s'/ //g'`"

if test $nproc = 4 ; then
	corecombs="1,2,4"
elif test $nproc = 6 ; then
	corecombs="1,2,4,6"
elif test $nproc = 8 ; then
	corecombs="1,2,4,8"
elif test $nproc = 12 ; then
	#corecombs="1,2,6,12"
        corecombs="1,2,4,6,8,12"
elif test $nproc = 16 ; then
	corecombs="1,2,4,8,16"
elif test $nproc = 24 ; then
	        corecombs="1,2,4,8"
elif test $nproc = 48 ; then
        # Istanbul fix
        corecombs="1,2,4,8,12"
elif test $nproc = 32 ; then
	# AIX fix
	corecombs="1,2,4,8,16"
elif test $nproc = 64 ; then
	# AIX fix
	corecombs="1,2,4,8,16"
else
	#corecombs=`seq $nproc| sed "s/ /,/g"`#evil : seq produces newlines
	# FIXME : does not work on weird SP5's
	corecombs=`seq $nproc| tr '\n' ' ' |  sed "s/\> \</,/g"`
#else
#	echo "uhm. did not recognize the number of cores!"
#	exit
fi
else
	corecombs=":"
fi


toytest=1 # NOTE: the only supported for now
RB=${RSBDB_RSBENCH}
#RB=./rsbench 
SW="--echo-arguments --verbose ${RSBDB_COMPARE_SWITCH} --times ${RSBDB_TIMES}"
RSB="-Fbo -qH -n $corecombs -R " 
RSBDB_SPMV="$RB -oa -Ob $SW $RSB" 
RSBDB_SPMVT="$RB -oa -Ob $SW $RSB --transpose" 
RSBDB_SPMVSYM="$RB -oa -Ob $SW $RSB --as-symmetric" 
RSBDB_SPSV="$RB -ot -Ob $SW $RSB"
OB_SPMV="$RB -oa -Ob $SW -Fo" 
CB_SPMV="$RB -oa -Ob $SW -Fb"

###############################################################################
if test x$toytest = x"1" ; then
###############################################################################
#seq=`seq 0 2`
seq=`seq 0 0`
spacings=`for p in $seq ; do echo $((2**p)) ; done  `
dseq="`seq $RSBDB_MINDIMLOG $RSBDB_MAXDIMLOG`"
#dseq=`seq 0 3`
#dims="512 1024"
dims=`for p in $dseq ; do if test $p = 0 ; then echo 1; else echo $((2**p)) $(( (3*2**p) / 2)) ; fi ; done  `
#incxs="1 2 4"
incxs="1"
###############################################################################
else
###############################################################################
#seq=`seq 14`
#seq=`seq 0 10`
seq=`seq 0 5`
spacings=`for p in $seq ; do echo $((2**p)) ; done  `
dims="512 1024 2048 4096"
incxs="1 2 4 8 16 32"
###############################################################################
fi
###############################################################################
incxs=${RSBDB_INCXS:-$incxs}
dims=${RSBDB_DIMS:-$dims}
spacings=${RSBDB_SPACINGS:-$spacings}
###############################################################################
# uncomment the following for testing purposes
#dims="256 512 "
#dims="$dims"
ddims="$dims"
sdims="$dims"
###############################################################################

tic() { echo `date +%Y%m%d%H%M%S`.`date +%s`; }

TYPES=`grep RSB_NUMERICAL_TYPE_PREPROCESSOR_SYMBOLS types.h | sed 's/  */ /g;s/"//g' | cut -d \  -f 3-`

for type in $TYPES; do  
mbln=microbench-librsb-$RSBDB_VERSION-$type
df=--dense
lf=--lower
# FIXME: AIX has problems recognizing long switches
# 20110109  but these short options seem buggy either way
#df=-d
#lf=-l

if test -z ${RSBDB_LABEL} ; then did=`tic` ; else did="${RSBDB_LABEL}" ; fi
export RSBDB_LISTFILE=${RSBDB_LISTFILE:-$mbln-$did-files.txt}
rm ${RSBDB_LISTFILE}

cp $0 $mbln-$did-script.sh
ldd ./rsbench  > $mbln-$did-ldd.txt
export RSBDB_FILESLIST="${RSBDB_FILESLIST} $mbln-$did-ldd.txt $mbln-$did-script.sh"

dfnopt=" -T $type --override-matrix-name "
if test -z ${RSBDB_LABEL} ; then did=`tic` ; else did="${RSBDB_LABEL}" ; fi
ln=$mbln-$did.tmp

for dim in $sdims ; do for spa in $spacings ; do 
	$RSBDB_SPMVSYM $lf $dim --generate-spacing $spa $dfnopt $dim
done ; done  2>&1 | tee $ln
if test -z ${RSBDB_LABEL} ; then dod=`tic` ; else dod="" ; fi
fn=${ln//.tmp/}$dod-symspmv-spaced.log
mv $ln $fn 
plotstuff "$fn"

if test -z ${RSBDB_LABEL} ; then did=`tic` ; else did="${RSBDB_LABEL}" ; fi
ln=$mbln-$did.tmp

for dim in $sdims ; do for spa in $spacings ; do 
	$RSBDB_SPMVT $df $dim --generate-spacing $spa $dfnopt $dim
done ; done  2>&1 | tee $ln
if test -z ${RSBDB_LABEL} ; then dod=`tic` ; else dod="" ; fi
fn=${ln//.tmp/}$dod-spmvt-spaced.log
mv $ln $fn 
plotstuff "$fn"

if test -z ${RSBDB_LABEL} ; then did=`tic` ; else did="${RSBDB_LABEL}" ; fi
ln=$mbln-$did.tmp

for dim in $sdims ; do for spa in $spacings ; do 
	$RSBDB_SPMV $df $dim --generate-spacing $spa $dfnopt $dim
done ; done  2>&1 | tee $ln
if test -z ${RSBDB_LABEL} ; then dod=`tic` ; else dod="" ; fi
fn=${ln//.tmp/}$dod-spmv-spaced.log
mv $ln $fn 
plotstuff "$fn"


if test -z ${RSBDB_LABEL} ; then did=`tic` ; else did="${RSBDB_LABEL}" ; fi
ln=$mbln-$did.tmp
for dim in $sdims ; do
	$RSBDB_SPSV $lf $dim $dfnopt $dim
done  2>&1 | tee $ln
if test -z ${RSBDB_LABEL} ; then dod=`tic` ; else dod="" ; fi
fn=${ln//.tmp/}$dod-spsv-solve.log
mv $ln $fn 
plotstuff "$fn"

if false ; then 

did=`tic`
ln=$mbln-$did.tmp
for dim in $sdims ; do
	$CB_SPMV $df $dim $dfnopt $dim
done  2>&1 | tee $ln
if test -z ${RSBDB_LABEL} ; then dod=`tic` ; else dod="" ; fi
fn=${ln//.tmp/}$dod-spmv-csr.log
mv $ln $fn 
plotstuff "$fn"


did=`tic`
ln=$mbln-$did.tmp
for dim in $sdims ; do
	$OB_SPMV $df $dim $dfnopt $dim
done  2>&1 | tee $ln
if test -z ${RSBDB_LABEL} ; then dod=`tic` ; else dod="" ; fi
fn=${ln//.tmp/}$dod-spmv-coo.log
mv $ln $fn 
plotstuff "$fn"


did=`tic`
ln=$mbln-$did.tmp
for dim in $ddims ; do for i in $incxs ; do
	$CB_SPMV $df $dim --incx $i $dfnopt $dim
done ; done  2>&1 | tee $ln
if test -z ${RSBDB_LABEL} ; then dod=`tic` ; else dod="" ; fi
fn=${ln//.tmp/}$dod-spmv-dense-csr-incx.log
mv $ln $fn 
plotstuff "$fn"


did=`tic`
ln=$mbln-$did.tmp
for dim in $ddims ; do for i in $incxs ; do
	$CB_SPMV $df $dim --incy $i $dfnopt $dim
done ; done  2>&1 | tee $ln
if test -z ${RSBDB_LABEL} ; then dod=`tic` ; else dod="" ; fi
fn=${ln//.tmp/}$dod-spmv-dense-csr-incy.log
mv $ln $fn 
plotstuff "$fn"

did=`tic`
ln=$mbln-$did.tmp
for dim in $ddims ; do for i in $incxs ; do
	$RSBDB_SPMV $df $dim --incx $i $dfnopt $dim
done ; done  2>&1 | tee $ln
if test -z ${RSBDB_LABEL} ; then dod=`tic` ; else dod="" ; fi
fn=${ln//.tmp/}$dod-spmv-dense-rsb-incx.log
mv $ln $fn 
plotstuff "$fn"


did=`tic`
ln=$mbln-$did.tmp
for dim in $ddims ; do for i in $incxs ; do
	$RSBDB_SPMV $df $dim --incy $i $dfnopt $dim
done ; done  2>&1 | tee $ln
if test -z ${RSBDB_LABEL} ; then dod=`tic` ; else dod="" ; fi
fn=${ln//.tmp/}$dod-spmv-dense-rsb-incy.log
mv $ln $fn 
plotstuff "$fn"

fi
done
#
echo ${RSBDB_FILESLIST} ${RSBDB_FILESLIST//log/gp} > ${RSBDB_LISTFILE}
#
# FIXME: missing with differing cache blocking, now

