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

# This script is horrible.
# This script is intended for the librsb developer usage.

#
# to force single thread execution non recursive, should invoke after exporting
# RSB_BENCH_FORMATS=o RSB_BENCH_THREADS=1 RSB_WANT_Q=" " RSB_BENCH_RECURSION=" "
#
# to force Z sorted COO:
# RSB_BENCH_FORMATS=o RSB_BENCH_THREADS=1 RSB_WANT_Q=" " RSB_BENCH_RECURSION=" " RSB_BENCH_EXTRA="--z-sorted-coo"
#
# to force Zig Zag CSR:
# RSB_BENCH_FORMATS=b RSB_BENCH_THREADS=1 RSB_WANT_Q=" " RSB_BENCH_RECURSION=" " RSB_BENCH_EXTRA="--zig-zag"
#
#svn revert INSTALL # so no svn diff
# this will fail often, but it's not a problem :)
#
# TODO: need to set affinity. shall work on:
# - OMP_CPU_BIND # 0 or 1
# - GOMP_CPU_AFFINITY='0-16:2' # alternate the CPUs
# - KMP_AFFINITY=verbose,granularity=thread,proclist=[0,1,2,3],explicit
# - KMP_AFFINITY=verbose,compact
# - KMP_AFFINITY=verbose,interleaved
# - KMP_AFFINITY=verbose,
# - KMP_VERSION=1
# 
# TODO: rename these variables
#
RSB_BENCH_THREADS=${RSB_BENCH_THREADS:=}
RSB_BENCH_HALFWORD=${RSB_WANT_Q:=}
RSB_BENCH_BB=${RSB_WANT_BB:=}
RSB_BENCH_EXTRA=${RSB_BENCH_EXTRA:=}
RSB_BENCH_RECURSION=${RSB_BENCH_RECURSION:=-R}
RSB_BENCH_FORMATS=${RSB_BENCH_FORMATS:=ob}
RSB_BENCH_COMPARE_OPTION=${RSB_BENCH_COMPARE_OPTION:=--compare-competitors}
#RSB_BENCH_AUTOTUNE_OPTION=${RSB_BENCH_AUTOTUNE_OPTION:=--want-autotune 3s10x-1t}
RSB_BENCH_AUTOTUNE_OPTION=${RSB_BENCH_AUTOTUNE_OPTION:=}
RSB_SKIP_SYMMETRIC=${RSB_SKIP_SYMMETRIC:=0}
RSB_SKIP_COMPLEX=${RSB_SKIP_COMPLEX:=0}
RSB_SKIP_UNSYMMETRIC=${RSB_SKIP_UNSYMMETRIC:=0}
RSB_MATDIR=${RSB_MATDIR:=}
RSB_TRANSPOSITIONS=${RSB_TRANSPOSITIONS:=--notranspose --transpose}
#RSB_TRANSPOSITIONS=${RSB_TRANSPOSITIONS:=--also-transpose}
RSB_BENCH_REPEAT_CONSTRUCTOR=${RSB_BENCH_REPEAT_CONSTRUCTOR:="5"}
RSB_WANT_CP_TO_SHM=${RSB_WANT_CP_TO_SHM:=0}
RSB_SHM=/dev/shm

#version=`svnversion` 
bid=`date +%Y%m%d%H%M%S`.`date +%s` # date id
version=`grep '\<VERSION\>' config.h | sed 's/.*VERSION *//g;s/"//g;s/0\.//g'`

# Matrices directory:
arg="$@"

#if true ; then
#if false ; then
if ! test -f rsbench ; then
#touch configure.ac
#sh autogen.sh
./configure --enable-sparse-blas-interface --enable-matrix-types=blas --enable-matrix-ops=blas --enable-openmp
make cleanall
make clean
make
make feedback || exit
[ -f rsbench ] || exit

fi

did=`date +%Y%m%d%H%M%S`.`date +%s` # date id
md=~/matrices

#matrices=$md/raefsky4.mtx
if true ; then
	matrices=""
	md=
	if test "x$arg" != x ; then
		md="$arg" ;echo $md;
	else
		if test x"${RSB_MATDIR}" != x ; then  md="${RSB_MATDIR}" ; fi
	fi
	om=`find $md -name \*.mtx -or -name \*.mtx.gz | sort`
	for m in $om ; do
		#f=$md/$m
		if test "x$RSB_SKIP_SYMMETRIC" = x1 ; then if head -n 1 "$m" | grep -i "symmetric\|hermitian" > /dev/null ; then echo "skipping $m" ; continue ; fi ; fi
		if test "x$RSB_SKIP_COMPLEX" = x1 ; then if head -n 1 "$m" | grep -i "complex" > /dev/null ; then echo "skipping $m" ; continue ; fi ; fi
		if test "x$RSB_SKIP_UNSYMMETRIC" = x1 ; then if head -n 1 "$m" | grep -i "general" > /dev/null ; then echo "skipping $m" ; continue ; fi ; fi
		f=$m
		matrices="$matrices $f"
	done
fi
#echo $matrices

hn=`hostname`
od=bench-svn$version-$hn-$did
mkdir -p $od || exit

# we gain processor info, if available.
cat /proc/cpuinfo > $od/$hn-cpuinfo.txt
cat /proc/meminfo > $od/$hn-meminfo.txt
x86info           > $od/$hn-x86info.txt
cpuid             > $od/$hn-cpuid.txt
cpuinfo           > $od/$hn-cpuinfo.txt
numactl --hardware > $od/$hn-numactl-H.txt # -H seems not to work on some systems
cp config.h         $od/ && gzip $od/config.h     
cp config.log       $od/ && gzip $od/config.log
cp Makefile         $od/ && gzip $od/Makefile
tar czf $od/sys-cpu.tgz /sys/devices/system/cpu/
cp bench/spmv.sh    $od/ && gzip $od/spmv.sh
env                > $od/env.txt

./rsbench -v  2>&1 > $od/rsbench-v.txt
./rsbench -C  2>&1 > $od/rsbench-C.txt
#./rsbench -M  2>&1 > $od/rsbench-M.txt # why comment this ? because it's very slow. (put it at the end)
./rsbench -I  2>&1 > $od/rsbench-I.txt
ldd ./rsbench  > $od/ldd-rsbench.txt

( for m in $matrices ; do echo $m ; done ) > $od/matrices-list.txt
( for m in $matrices ; do ls -l $m ; done ) > $od/matrices-list-l.txt

#nproc="`cat /proc/cpuinfo | grep ^processor | wc -l`"
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

corecombs=${RSB_BENCH_THREADS:=$corecombs}


if test x`hostname` = x"drachma" ; then n="numactl --physcpubind=0-11 --membind=0-1 " ; else n="" ; fi

for take in '-take1' '-take2' ;
do
ont=$nproc
#ont=""

#export OMP_NUM_THREADS=$ont
#types=${RSB_BENCH_TYPES:=D S C Z}
types=${RSB_BENCH_TYPES:=D}
if test "x$RSB_CACHE_FLUSH" = x1 ; then cache_flush=--cache-flush ; else cache_flush=''; fi

#for s in   -R ; do for F in br bc ; do for m in $matrices ; do
#for s in  -R ; do for F in br ; do for m in $matrices ; do
for T in $types ; do for s in "$RSB_BENCH_RECURSION" ; do for F in "$RSB_BENCH_FORMATS" ; do for m in $matrices ; do 
	mn=`basename $m`
	if test x${RSB_WANT_CP_TO_SHM} = x"1" ; then
		cp --dereference $m ${RSB_SHM} || exit ; m="${RSB_SHM}/$mn" ;
	fi

	for transa in ${RSB_TRANSPOSITIONS} ; do

#for s in   -R ; do for F in br bc ; do for m in $md/*.mtx ; do
	times=100
	# IBM sp5's stat is BUGGED : do not use it
        if test $(( `du -sk "$m" | cut -f 1` > 500000 )) = 1 ; then times=$((times/2)) ; fi
#	if test $(( `du -sb "$m" | cut -f 1` > 500000000 )) = 1 ; then times=$((times/2)) ; fi
#	if test $(( `stat --format "%s" "$m"` > 500000000 )) = 1 ; then times=$((times/2)) ; fi
	did=`date +%Y%m%d%H%M%S`.`date +%s` # date id
	# FIXME : using -TD segfaults (it's due to the 'T' switch)!
	# use short options, as getopt_long does not exist on some platforms.

	DQ="-qH" # --cache-flush
	Q=${RSB_BENCH_HALFWORD:-$DQ}
	BB=${RSB_BENCH_BB:- --bounded-box=1}
	Q="$Q $BB -V $transa --repeat-constructor $RSB_BENCH_REPEAT_CONSTRUCTOR " # --cache-flush

	cmd="./rsbench -oa -Ob -F$F -f $m $s -t $times -r1 -c1 ${RSB_BENCH_COMPARE_OPTION} $Q -n $corecombs -T $T $cache_flush $RSB_BENCH_EXTRA"
	echo $cmd
	$cmd 2>&1 | tee $od/bench-spmv-svn$version-$hn-$did-$mn-$s$F-T$T-omp-"$ont"''$transa'-'cores$take.log

done ;
	if test x${RSB_WANT_CP_TO_SHM} = x"1" ; then
		rm $m || exit
	fi
done ; done ; done ; done

done

if [ `hostname | cut -c 1-6` != "helios" ]; then # 20120413 temporary!
./rsbench -M  2>&1 > $od/rsbench-M.txt # why comment this ? because it's very slow. (put it at the end)
fi
./rsbench -F  2>&1 > $od/rsbench-F.txt # 
did=`date +%Y%m%d%H%M%S`.`date +%s` # date id
echo "benchmarking began at $bid"
echo "benchmarking ended at $did"
mv $od $od-$did || exit



