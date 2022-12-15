#!/bin/bash
#
# Copyright (C) 2008-2017 Michele Martone
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
# This script configures and builds binaries for a number of different setups.
# It is intended to work on the developer machine.

# TO DO:
# Need cryptographic hashes.
#
RSBB_DIR=${RSBB_DIR:-$HOME/src/src-tmp/} # an absolute address, please
RSBB_URL=${RSBB_URL:-svn+ssh://user@host/repository/trunk}
#RSBB_SVN=echo
#RSBB_SVN=svn
RSBB_MAKE=${RSBB_MAKE:-make}
RSBB_SVN=${RSBB_SVN:-svn}
RSBB_SVN_OPTIONS=${RSBB_SVN_OPTIONS:-}
RSBB_RUN_CONFIGURE=${RSBB_RUN_CONFIGURE:-1}
RSBB_RUN_AUTOGEN=${RSBB_RUN_AUTOGEN:-1}
RSBB_RUN_MAKE=${RSBB_RUN_MAKE:-1}
RSBB_RUN_TOUCH_HEADERS_TO_OLD=${RSBB_RUN_TOUCH_HEADERS_TO_OLD:-0}
RSBB_RUN_MAKE_CLEAN=${RSBB_RUN_MAKE_CLEAN:-1}
RSBB_RUN_MAKE_CLEANALL=${RSBB_RUN_MAKE_CLEANALL:-0}
RSBB_RUN_MAKE_TESTRUN=${RSBB_RUN_MAKE_TESTRUN:-1}
RSBB_RUN_MAKE_TESTS=${RSBB_RUN_MAKE_TESTS:-0}
RSBB_EXIT_ON_FAILED_TESTS=${RSBB_EXIT_ON_FAILED_TESTS:-1}
RSBB_RUN_MAKE_QTESTS=${RSBB_RUN_MAKE_QTESTS:-1}
RSBB_RUN_MAKE_QQTESTS=${RSBB_RUN_MAKE_QQTESTS:-0}
RSBB_RUN_MAKE_JOBS=${RSBB_RUN_MAKE_JOBS:-2}
RSBB_RUN_BUILD_ARCHIVE=${RSBB_RUN_BUILD_ARCHIVE:-1}
RSBB_RUN_MAKE_DOX=${RSBB_RUN_MAKE_DOX:-0}
RSBB_RUN_MAKE_DOXONLY=${RSBB_RUN_MAKE_DOXONLY:-0}
RSBB_RUN_MAKE_INSTALL=${RSBB_RUN_MAKE_INSTALL:-0}
RSBB_RUN_MAKE_LOCAL_LIBRSB_LINK=${RSBB_RUN_MAKE_LOCAL_LIBRSB_LINK:-0}
RSBB_WANT_LIST_FILE=${RSBB_WANT_LIST_FILE:-1}
RSBB_WANT_README_FILE=${RSBB_WANT_README_FILE:-1}
RSBB_REBUILD_STATIC_RSBENCH=${RSBB_REBUILD_STATIC_RSBENCH:-0}
RSBB_RUN_SVN_CO=${RSBB_RUN_SVN_CO:-1}
RSBB_RUN_SVN_UP=${RSBB_RUN_SVN_UP:-1}
RSBB_RUN_SVN_DIFF=${RSBB_RUN_SVN_DIFF:-0}
RSBB_RUN_SVN_REVERT=${RSBB_RUN_SVN_REVERT:-0}
RSBB_RUN_SVN_CLEANUP=${RSBB_RUN_SVN_CLEANUP:-0}
RSBB_DIST_URL=${RSBB_DIST_URL:-/dev/null}
RSBB_DIST_DOWNLOAD=${RSBB_DIST_DOWNLOAD:-0}
RSBB_DIST_UNPACK=${RSBB_DIST_UNPACK:-0}
RSBB_DIST_UNPACK_ARCHIVE=${RSBB_DIST_UNPACK_ARCHIVE:-$RSBB_DIR/librsb.tar.gz}
RSBB_BASENAME=${RSBB_BASENAME:-librsb-batch-static}
RSBB_INSTALL_PREFIX=${RSBB_INSTALL_PREFIX:-/usr/local/}
RSBB_WANT_OVERRIDE_PREFIX_WITH_SAME_DIR=${RSBB_WANT_OVERRIDE_PREFIX_WITH_SAME_DIR:-0}
RSBB_WANT_SKIP_CANONICAL_CFLAGS=${RSBB_WANT_SKIP_CANONICAL_CFLAGS:-0}
RSBB_CONFIGURE_ALTERNATIVES=${RSBB_CONFIGURE_ALTERNATIVES:---disable-openmp --enable-openmp}
RSBB_DIST_UNPACK_TARTRANSFEXP=${RSBB_DIST_UNPACK_TARTRANSFEXP:-'s/librsb-*[^\/]*\///'} # was: 's/librsb-*[0-9:M.]\+\///' 
RSBB_CONFIGURE_ADD=${RSBB_CONFIGURE_ADD:-}
RSBB_RUN_PRECONFIGURE_HOOK=${RSBB_RUN_PRECONFIGURE_HOOK:-}
RSBB_RUN_PREMAKE_HOOK=${RSBB_RUN_PREMAKE_HOOK:-}
RSBB_RUN_POSTMAKE_HOOK=${RSBB_RUN_POSTMAKE_HOOK:-}
RSBB_MAKE_EXTRA_CFLAGS=${RSBB_MAKE_EXTRA_CFLAGS:-}
CC=${CC:-gcc}
FC=${FC:-gfortran}
CXX=${CXX:-g++}
RSBB_CC_ALTERNATIVES=${RSBB_CC_ALTERNATIVES:-$CC}
RSBB_CFLAGS_ALTERNATIVES=${RSBB_CFLAGS_ALTERNATIVES:-}
RSBB_CFLAGS_ADD=${RSBB_CFLAGS_ADD:--fPIC}
RSBB_MKL=${RSBB_MKL:-}
RSBB_WANT_TYPES=${RSBB_WANT_TYPES:-}
RSBB_WANT_TOLERATE_MAKE_FAILURE=${RSBB_WANT_TOLERATE_MAKE_FAILURE:-0}
#RSBB_WANT_TYPES=${RSBB_WANT_TYPES:---enable-matrix-types=double}
#RSBB_WANT_TYPES=${RSBB_WANT_TYPES:---enable-matrix-types=double,float,float complex,double complex}
#RSBB_WANT_TYPES=${RSBB_WANT_TYPES//%/ }
RSBB_DIR_NAME_APPEND=${RSBB_DIR_NAME_APPEND:-}
RSBB_DIR_NAME_PRETAG=${RSBB_DIR_NAME_PRETAG:-}
RSBB_DIR_NAME_PREPEND=${RSBB_DIR_NAME_PREPEND:-}
DATE=date
#
if test x$# = x0 ; then 
	echo "Usage: $0 ''"
	echo "Any of the following variables may be set:"
	echo
	declare | grep RSBB_
	echo
	exit
else
	true ;
fi
ECHO=echo
ECHO=
CP=cp
BBL=batchbuild.log
#cct='--enable-matrix-types=double,double complex'
#cct=--disable-openmp
cct=
ARCH=`uname -m`
RSBB_BN=${RSBB_BASENAME}-${ARCH}
LIST=${RSBB_DIR}/${RSBB_BN}-files.txt
RSBB_RF=
if test x$RSBB_WANT_README_FILE = x1 ; then
	RSBB_RF=${RSBB_DIR}/${RSBB_BN}-README.txt
fi
BA=${RSBB_DIR}/${RSBB_BN}.tar.gz
sf="-all-static"
#for co in "--enable-openmp"
#for co in "--enable-openmp" "--disable-openmp"
rm -f ${LIST}
#
mkdir -p $RSBB_DIR || exit -1
#
cd $RSBB_DIR || exit -1
#
if test x$RSBB_WANT_README_FILE = x1 ; then
	echo "README (log) file for librsb build." > $RSBB_RF
	#echo "$USER@$HOSTNAME" >> $RSBB_RF
	echo "Host: $HOSTNAME" >> $RSBB_RF
	uname -a >> $RSBB_RF
	date >> $RSBB_RF
	echo "Relevant built variables: " >> $RSBB_RF
	declare | grep RSBB_  >> $RSBB_RF
	echo " " >> $RSBB_RF
fi
#
if test x${RSBB_DIST_DOWNLOAD} = x1 ; then
	$ECHO wget ${RSBB_DIST_URL} -O librsb.tar.gz || exit -1
fi
#
#for co in "--enable-openmp" "--disable-openmp"
for cc in $RSBB_CC_ALTERNATIVES
do
export CC="${cc}"
#if test x$cc = xicc ; then export FC="ifort"; fi
if test x`basename $cc` = xicc ; then export FC=`dirname $cc`/"ifort"; fi
for co in $RSBB_CONFIGURE_ALTERNATIVES
#for co in  "--disable-openmp"
do
alternative_already_done=no # when RSBB_CFLAGS_ALTERNATIVES is one of the presets
for bo in "-O0 -ggdb" "-O2 -pg" "-O3" "${RSBB_CFLAGS_ALTERNATIVES}"
#for bo in "-O3"  "${RSBB_CFLAGS_ALTERNATIVES}"
#for bo in "-O0 -ggdb"  "${RSBB_CFLAGS_ALTERNATIVES}"
#for bo in "-O2 -pg" "-O3"
#for bo in "-O0 -ggdb"
do
	if test x"$bo" = x"$RSBB_CFLAGS_ALTERNATIVES" -a x"$bo" = x ; then continue ; fi
	if test x"1" = x"$RSBB_WANT_SKIP_CANONICAL_CFLAGS" -a x"$bo" != x"$RSBB_CFLAGS_ALTERNATIVES" ; then continue ; fi
	if test x"$alternative_already_done" = x"yes" ; then continue; fi
	if test x"$bo" = x"$RSBB_CFLAGS_ALTERNATIVES" ; then alternative_already_done=yes ; fi
	#bo="$bo ${RSBB_MAKE_EXTRA_CFLAGS}"
	$DATE
	#TAG=${bo// /}-${co// /}
	TAG=`basename ${CC}`-${bo// /}-${co// /}
	BT=${RSBB_DIR_NAME_PREPEND}$RSBB_BN-${RSBB_DIR_NAME_PRETAG}${TAG}${RSBB_DIR_NAME_APPEND}
	BD=$RSBB_DIR/${BT}
	CFLAGS="${bo} ${RSBB_CFLAGS_ADD} ${RSBB_MAKE_EXTRA_CFLAGS}"
	export CFLAGS="$CFLAGS"
	export CXXFLAGS=${CXXFLAGS:-$CFLAGS}
	#export CXXFLAGS="$CFLAGS"
	export CXX
	# $ECHO $RSBB_SVN info $RSBB_URL || exit -1
	#
	if test "x${RSBB_DIST_UNPACK}" = x1 ; then
		mkdir -p librsb/ || exit -1 
		tar --transform "${RSBB_DIST_UNPACK_TARTRANSFEXP}" -C librsb -xzf ${RSBB_DIST_UNPACK_ARCHIVE} || exit -1
		#cp -fvRp librsb ${BD}|| exit -1 
		#mv librsb ${BD} || exit -1 
		mkdir -p ${BD} || exit -1 
		#mv -f librsb/* ${BD}/ || exit -1 
		rsync -avz librsb/ ${BD}/ || exit -1 
		rm -fR librsb  || exit -1 
		#rmdir librsb  || exit -1 
	fi
	#
	if test "x${RSBB_RUN_SVN_CO}" = x1 ; then
		$ECHO $RSBB_SVN ${RSBB_SVN_OPTIONS} --force co ${RSBB_URL} ${BD} || exit -1
	fi
	#
	if test "x${RSBB_RUN_SVN_DIFF}" = x1 ; then
		$ECHO $RSBB_SVN ${RSBB_SVN_OPTIONS} info ${RSBB_URL} || exit -1
		$ECHO $RSBB_SVN ${RSBB_SVN_OPTIONS}         diff $BD || exit -1
	fi
	#
	if test "x${RSBB_RUN_SVN_REVERT}" = x1 ; then
		$ECHO $RSBB_SVN ${RSBB_SVN_OPTIONS} info ${RSBB_URL} || exit -1
		$ECHO $RSBB_SVN ${RSBB_SVN_OPTIONS}         revert `$RSBB_SVN ${RSBB_SVN_OPTIONS}         ls $BD` || exit -1
	fi
	#
	if test "x${RSBB_RUN_SVN_CLEANUP}" = x1 ; then
		$ECHO $RSBB_SVN ${RSBB_SVN_OPTIONS} info ${RSBB_URL} || exit -1
		$ECHO $RSBB_SVN ${RSBB_SVN_OPTIONS}         cleanup $BD || exit -1
	fi
	#
	if test "x${RSBB_RUN_SVN_UP}" = x1 ; then
		$ECHO $RSBB_SVN ${RSBB_SVN_OPTIONS} info ${RSBB_URL} || exit -1
		$ECHO $RSBB_SVN ${RSBB_SVN_OPTIONS} --force up  $BD || exit -1
	fi
	#
	echo will cd $BD || exit -1
	$ECHO cd $BD || exit -1
	{
	if test x"$RSBB_RUN_PRE_HOOK" != x ; then
		eval "$RSBB_RUN_PRE_HOOK"
	fi
	#if test x"$RSBB_RUN_MAKE" = x1 -a '!' -e Makefile ; then
	if test '!' -e Makefile && set | grep 'RSBB_RUN_MAKE.*=1$' ; then
		RSBB_RUN_CONFIGURE=1 # create a Makefile if it does not exist
	fi
	if test x"$RSBB_RUN_CONFIGURE" = x1 ; then
		if test configure.ac -nt configure ; then $ECHO sh autogen.sh || exit -1 ; fi
	fi
	mkdir -p m4 # for Makefile.am
	if test x"$RSBB_RUN_AUTOGEN" = x1 ; then
		sh autogen.sh || exit -1
		# FIXME: temporary, for Helios:
		#libtoolize -c # || exit -1
		autoreconf --install --force || true # 20121016
	fi
	if test x"$RSBB_WANT_OVERRIDE_PREFIX_WITH_SAME_DIR" = x1 ; then
		RSBB_INSTALL_PREFIX="`pwd`/local/"
	fi
	cco="--prefix=${RSBB_INSTALL_PREFIX}"
	#cco='--prefix=/usr/local/'
#
	if test x"$RSBB_RUN_PRECONFIGURE_HOOK" != x ; then
		# FIXME; seems like the following export does not work
		#export RSBB_CONFIGURE_ADD="${RSBB_CONFIGURE_ADD} $co" # seems incorrect
		coo="$co ${RSBB_CONFIGURE_ADD}" # seems correct
		#export co="$co ${RSBB_CONFIGURE_ADD}" # seems correct
		#echo "RSBB_CONFIGURE_ADD is now: ${RSBB_CONFIGURE_ADD}"
		export CFLAGS="${CFLAGS}";
		eval "$RSBB_RUN_PRECONFIGURE_HOOK"
		if test x"$RSBB_RUN_CONFIGURE" = x1 ; then
			echo "configure will now be invoked with: ${coo}"
		fi
	else
		coo="$co " # seems correct
	fi
#
	RSBB_MKL_OPTION=
	if test x"$RSBB_MKL" != x ; then
		RSBB_MKL_OPTION=--with-mkl="$RSBB_MKL"
	fi
	if test x"$RSBB_RUN_MAKE_LOCAL_LIBRSB_LINK" = x1 ; then
		$ECHO rm -f local  || exit -1
		$ECHO ln -v -f -s `echo ${BD} | sed s/${RSBB_BASENAME}/librsb/g`/local local  || exit -1
	fi
	if test x"$RSBB_RUN_CONFIGURE" = x1 ; then
		$ECHO ./configure $coo "${RSBB_WANT_TYPES}" "${RSBB_MKL_OPTION}" CFLAGS="${CFLAGS}" ${cco} ${cct} CC=${CC} FC=${FC} CXX=${CXX} \
			LIBRSB_CONFIG=`pwd`/local/bin/librsb-config # this is an extra for sparsersb
	fi
	if test x"$RSBB_RUN_TOUCH_HEADERS_TO_OLD" = x1 ; then
		# FIXME: is.h is not necessary old :P
		#touch config.h rsb.h -r is.h || exit -1
		#touch *.h *.m4 -r is.h || exit -1
		touch *.h -r is.h || exit -1
	fi
	if test x"$RSBB_RUN_MAKE_CLEAN" = x1 ; then
		$ECHO ${RSBB_MAKE} clean || exit -1
	fi
	if test x"$RSBB_RUN_MAKE_CLEANALL" = x1 ; then
		$ECHO ${RSBB_MAKE} cleanall || exit -1
	fi
	if test x"$RSBB_RUN_MAKE" = x1 ; then
		if test x"$RSBB_RUN_PREMAKE_HOOK" != x ; then
			eval "$RSBB_RUN_PREMAKE_HOOK"
		fi
		RSBB_MAKE_ERROR_EXIT=exit
		if test x"$RSBB_WANT_TOLERATE_MAKE_FAILURE" = x1 ; then
			RSBB_MAKE_ERROR_EXIT=echo
		fi
		$ECHO ${RSBB_MAKE} all -j ${RSBB_RUN_MAKE_JOBS} # || exit -1
		#$ECHO ${RSBB_MAKE} all -j 1 || ${RSBB_MAKE_ERROR_EXIT} -1 # FIXME: this is temporary
		if test x"$RSBB_RUN_POSTMAKE_HOOK" != x ; then
			eval "$RSBB_RUN_POSTMAKE_HOOK"
		fi
	fi
	if test x"$RSBB_RUN_MAKE_QTESTS" = x1 ; then
		$ECHO ${RSBB_MAKE} qtests || exit -1
	fi
	if test x"$RSBB_RUN_MAKE_QQTESTS" = x1 ; then
		$ECHO ${RSBB_MAKE} qqtests || exit -1
	fi
	if test x"$RSBB_RUN_MAKE_TESTS" = x1 ; then
		$ECHO ${RSBB_MAKE} tests || if test x${RSBB_EXIT_ON_FAILED_TESTS} != x0 ; then exit -1 ; fi # too slow
	fi
	if test x"$RSBB_RUN_MAKE_DOX" = x1 ; then
		$ECHO ${RSBB_MAKE} dox || exit -1
	fi
	if test x"$RSBB_RUN_MAKE_DOXONLY" = x1 ; then
		$ECHO ${RSBB_MAKE} doxonly || exit -1
	fi
	if test x"$RSBB_RUN_MAKE_INSTALL" = x1 ; then
		$ECHO ${RSBB_MAKE} install || exit -1
	fi
	if test x"$RSBB_REBUILD_STATIC_RSBENCH" = x1 ; then
	if test x"$RSBB_RUN_MAKE" = x1 ; then
		echo "Will attempt to rebuild statically.."
		$ECHO rm -f rsbench 
		#$ECHO ${RSBB_MAKE} CFLAGS="${CFLAGS} ${sf}" || exit -1
		$ECHO ${RSBB_MAKE} AM_CFLAGS="${CFLAGS} ${sf}" || exit -1
	fi
	fi
	if test x"$RSBB_RUN_MAKE_TESTRUN" = x1 ; then
		./rsbench || exit -1
	fi
	if test x"$RSBB_WANT_LIST_FILE" = x1 ; then
		sed 's/^prefix=.*$/prefix="\/usr\/local\/"/g' `pwd`/librsb-config > `pwd`/librsb-config-local
		#ls -ltr 	`pwd`/librsb-config
		#ls -ltr 	`pwd`/librsb-config-local
		#for FN in rsbench librsb.a rsb.h rsb-spblas.h librsb-config rsb-types.h ; do
		for FN in rsbench librsb.a rsb.h rsb-spblas.h librsb-config-local rsb-types.h ; do
			#echo "$BD/$FN" >> ${LIST}
			echo "$BT/$FN" >> ${LIST}
			#echo "$BT/local/$FN" >> ${LIST}
		done || exit -1
	fi
	#$ECHO $CP rsbench rsbench-$TAG || exit -1
	#$ECHO $CP librsbench.a librsbench-$TAG.a || exit -1
	#$ECHO $CP librsb-config librsb-config-$TAG || exit -1
	} # > $BD/$BBL
	# TODO: place email report sending here, or later on
	$DATE
	# 
	$ECHO cd -
	if test x"$RSBB_RUN_POST_HOOK" != x ; then
		eval "$RSBB_RUN_POST_HOOK"
	fi
done
done
done
	cd $RSBB_DIR || exit -1
	if test x$RSBB_WANT_LIST_FILE = x1 -a x$RSBB_RUN_BUILD_ARCHIVE = x1 ; then
		md5sum `cat ${LIST}` > ${LIST}.md5
		tar cvzf ${BA} --transform='s/-local//' `cat ${LIST}` # ${LIST}.md5
		#tar cvzf ${BA} `cat ${LIST}` # ${LIST}.md5
		md5sum ${BA} > ${BA}.md5
		echo "contents:"
		tar tzvf  ${BA}
		ls -ltr  ${BA}
	fi
	if test x$RSBB_WANT_README_FILE = x1 ; then
		echo "Build terminated at: " >> $RSBB_RF
		date >> $RSBB_RF
	fi
exit
