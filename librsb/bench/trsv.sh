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

#svn revert INSTALL # so no svn diff
# this will fail often, but it's not a problem :)

#version=`svnversion` 
version=`grep '\<VERSION\>' config.h | sed 's/.*VERSION *//g;s/"//g;s/0\.//g'`

arg="$1"

#if true ; then
#if false ; then
if ! test -f rsbench ; then
#touch configure.ac
#sh autogen.sh
./configure --enable-sparse-blas-interface --enable-matrix-types=blas --enable-matrix-ops=blas --disable-vb --with-openmp
make cleanall
make clean
make
make feedback || exit
[ -f rsbench ] || exit

fi

did=`date +%Y%m%d%H%M%S`.`date +%s` # date id
md=~/matrices_csb
md=~/matrices_hawaii
md=~/matrices

matrices=$md/raefsky4.mtx
# we use L factors of 
if false ; then
	md=~/matrices
	om="memplus.mtx wang4.mtx ex11.mtx raefsky4.mtx goodwin.mtx lhr10.mtx"
	for m in $om ; do
		f=$md/$m
		echo $f '->' ${f/.mtx/-superlu-L.mtx}
		./rsbench -L $f > ${f/.mtx/-superlu.mtx}
	done
else
	matrices=""
	md=~/matrices
	[ -d ~/matrices_ussv_colamd_L ] && md=~/matrices_ussv_colamd_L
	[ -d /scratch/downloads/matrices_ussv_colamd_L/ ] && md=/scratch/downloads/matrices_ussv_colamd_L/
	[ -d /opt/martone/matrices_ussv_colamd_L/ ] && md=/opt/martone/matrices_ussv_colamd_L/
	if test "x$arg" != x ; then md="$arg" ;echo $md; fi
	

	om=`find $md -name \*.mtx`
	#om="memplus.mtx wang4.mtx ex11.mtx raefsky4.mtx goodwin.mtx lhr10.mtx"
	#mom="av41092 FEM_3D g7jac180 g7jac200 garon2 ns3Da ohne2 para-9 poisson3Db rajat31 rma10 sme3Dc torso1 venkat50"
	#xom="stomach twotone"
	for m in $om ; do
		#f=$md/$m
		f=$m
		#ls -l $f
#		echo $f
	#	echo $f '->' ${f/.mtx/-superlu-L.mtx}
	#	matrices="$matrices "${f/.mtx/-superlu.mtx}
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
cp bench/trsv.sh    $od/ && gzip $od/trsv.sh
env                > $od/env.txt

./rsbench -v  2>&1 > $od/rsbench-v.txt
./rsbench -C  2>&1 > $od/rsbench-C.txt
./rsbench -M  2>&1 > $od/rsbench-M.txt
./rsbench -F  2>&1 > $od/rsbench-F.txt
./rsbench -I  2>&1 > $od/rsbench-I.txt
ldd ./rsbench  > $od/ldd-rsbench.txt

#nproc="`cat /proc/cpuinfo | grep ^processor | wc -l`"
nproc="`./rsbench -I | grep processors.online | sed 's/.*://'`"

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
elif test $nproc = 32 ; then
	# AIX fix
	corecombs="1,2,4,8,16"
elif test $nproc = 48 ; then
        # Istanbul fix
        corecombs="1,2,4,8,12"
elif test $nproc = 64 ; then
	# AIX fix
	corecombs="1,2,4,8,16"
else
	#corecombs=`seq $nproc| sed "s/ /,/g"`#evil : seq produces newlines
	# FIXME : does not work on weird SP5's (nproc keeps having a heading,useless space, and therefore we have to quote later on)
	corecombs=`seq $nproc| tr '\n' ' ' |  sed "s/\> \</,/g"`
fi

for take in '-take1' '-take2' ;
do
#for ont in `seq $nproc` ;
#do
nproc=${nproc/ /}
ont=${nproc/ /}
 
repeat_constructor=${repeat_constructor:="1"}

if test x`hostname` = x"drachma" ; then n="numactl --physcpubind=0-11 --membind=0-1 " ; else n="" ; fi

#export OMP_NUM_THREADS=$ont

( for m in $matrices ; do echo $m ; done ) > $od/matrices-list.txt
( for m in $matrices ; do ls -l $m ; done ) > $od/matrices-list-l.txt

#export OMP_NUM_THREADS=$ont
#types=${RSB_BENCH_TYPES:=D S C Z}
types=${RSB_BENCH_TYPES:=D}
if test "x$RSB_CACHE_FLUSH" = x1 ; then cache_flush=--cache-flush ; else cache_flush=''; fi

#for s in   -R ; do for F in br bc ; do for m in $matrices ; do
for T in $types ; do for s in   -R ; do for F in ob ; do for m in $matrices ; do for transa in --transpose --notranspose ; do
#for s in   -R ; do for F in br bc ; do for m in $md/*.mtx ; do
	mn=`basename $m`
	times=100
	# IBM sp5's stat is BUGGED : do not use it
#	if test $(( `du -sb "$m" | cut -f 1` > 500000000 )) = 1 ; then times=$((times/2)) ; fi
	# AIX does not have -b. but it has -k indeed.
	if test $(( `du -sk "$m" | cut -f 1` > 500000 )) = 1 ; then times=$((times/2)) ; fi
#	if test $(( `stat --format "%s" "$m"` > 500000000 )) = 1 ; then times=$((times/2)) ; fi
	did=`date +%Y%m%d%H%M%S`.`date +%s` # date id

	DQ="-qH" # --cache-flush
	Q=${RSB_WANT_Q:-$DQ}
	BB=${RSB_WANT_BB:- --bounded-box=1}
	Q="$Q $BB -V $transa --repeat-constructor $repeat_constructor" # --cache-flush

	# FIXME : using -TD segfaults (it's due to the 'T' switch)!
	# use short options, as getopt_long does not exist on some platforms.
	#cmd="$n ./rsbench -ot -Ob -F$F -f $m $s -t $times -r1 -c1 -Q -n $corecombs"
	cmd="$n ./rsbench -ot -Ob -F$F -f $m $s -t $times -r1 -c1 --compare-competitors $Q -n $corecombs -T $T $cache_flush --only-lower-triangle"
	#cmd="./rsbench -ot -Ob -F$F -f $m $s -t $times -r1 -c1 -qCDL -n $corecombs"
	echo $cmd
	$cmd 2>&1 | tee "$od/bench-trsv-svn$version-$hn-$did-$mn-$s$F-T$T-omp-"$ont"''$transa'-'''cores$take.log"

done ; done ; done ; done ; done

done
#done
did=`date +%Y%m%d%H%M%S`.`date +%s` # date id
mv $od $od-$did || exit



