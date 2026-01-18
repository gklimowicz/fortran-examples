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

# This script is intended for the librsb developer usage.

set -e
set -x
if grep RSB_REINIT_SINGLE_VALUE *.c --exclude rsb_rsb.c ; then exit 255 ; fi
#if grep '//' *.c  ; then exit 255 ; fi # TODO: activate this.
if grep -n 'RSB_DO_ERR_RETURN\>' rsb_rsb.c  ; then exit 255; else true ; fi
if cpp rsb.h | grep '()$'   ; then echo '[!] failed'; exit 255 ; else true ; fi
#for f in *.h ; do if cpp $f | grep '()$'   ; then echo '[!] failed'; exit 255 ; else true ; fi ; done
if grep -n --exclude=rsb_rsb.c 'RSB_DO_ERR_RETURN_INTERFACE\>' *.c ; then exit 255 ; else true ; fi
for RSBLIB in .libs/*.a; do
	test -f ${RSBLIB}
	if nm ${RSBLIB}  | grep '\s[DG]\s' | grep -v '\s[DG]\s''rsb_' | grep -v '\sD\s_' ; then exit 255 ; else true ; fi
	if nm ${RSBLIB}  | grep '\<T\>' | sed 's/^.*\s//g' | grep -v '^\(rsb\|BLAS\|blas\|__\)' ; then exit 255 ; else true ; fi
	if nm ${RSBLIB}  | grep ' C [a-z]' | grep -v ' C rsb_' ; then exit 255 ; else true ; fi
	if test `nm ${RSBLIB} | grep -v '\<T\> rsb__'  | grep '\<T\> rsb_' | cut -d \  -f 3 | wc -l` -gt 173 ; then exit 255 ; else true; fi
done
if ar t .libs/librsb.a | grep -v  '_l\?a-rsb\|^rsbpp\|^rsb\.o$' | grep -v ^rsb_ ; then echo '[!] failed source filenames check'; exit 255; else true ; fi
test -n "${abs_top_srcdir}"
cd ${abs_top_srcdir}
if which flawfinder; then
	test -f rsb_rsb.c
	flawfinder rsb_rsb.c | tee flawfinder.log
	echo "output of running flawfinder in rats.log"
fi
if which rats; then
	test -f rsb_rsb.c
	rats rsb_rsb.c | tee rats.log
	echo "output of running rats in rats.log"
fi
if cat *.F90 examples/*.F90| sed 's/!.*$//g'| grep '^.\{73,\}' ; then echo 'Some source code exceeds 72 chars!'; fi
if grep 'x"ac_cv_' configure.ac */configure.ac; then false ; else true; fi
if cat examples/*.cpp| grep '^.\{155,\}' ; then echo 'Some source code reaches 155 chars!'; fi
test -f NEWS
test -f README
cd - # back into ${abs_top_builddir}
( ! find * | grep "home.*$USER" )
# ( ! find * | grep "$HOSTNAME" )
# ( ! find doc/* -exec grep -l -i "$HOSTNAME" '{}' ';' -print; )
if grep -n '[^ ]\\leftarrow' doc/Doxyfile ; then exit 255; else true ; fi
if grep -n '^[^	].\{80,\}' README ; then exit 255 ; else true ; fi
#if grep -n '^.\{81,\}' NEWS       ; then exit 255 ; else true ; fi
if grep '	' *.F90 */*.F90 ; then exit 255 ; else true ; fi
for f in  examples/*0 examples/*c   ;
	do grep -q  -F  $f rsb.h || echo "$f is missing from among rsb.h mentioned examples!";
done

OUT=`./librsb-config --I_opts`
test `echo "$OUT" |  wc -l` -le 1
OUT=`./librsb-config --cxxflags`
test `echo "$OUT" |  wc -l` -le 1
OUT=`./librsb-config --ldflags`
test `echo "$OUT" |  wc -l` -le 1
OUT=`./librsb-config --libs --extra_libs`
test `echo "$OUT" |  wc -l` -le 1

./rsbench -E 0.1s || exit 255
#
if test -z "${RSB_DEVTESTS_NO_LIMITS_TEST}" ; then
	if nm ./rsbench | grep asan_init ; then
		# https://github.com/google/sanitizers/wiki/SanitizerCommonFlags
		export LSAN_OPTIONS=verbosity=1:allocator_may_return_null=true                 # https://github.com/google/sanitizers/wiki/ThreadSanitizerFlags
		export ASAN_OPTIONS=verbosity=1:allocator_may_return_null=true:halt_on_error=0 # https://github.com/google/sanitizers/wiki/AddressSanitizerFlags
	fi
	set +e
	./rsbench --limits-testing
	RC="$?"
	set -e
	test -n "$RC" # nonempty return code?
	if test "$RC" = 0; then
		true
	elif test "$RC" -gt 128; then
		echo "#Limits testing was killed by signal $(($RC-128)) (SIG`kill -l $RC`))!"
		if nm ./rsbench | grep asan_init ; then
			echo "Seems like -lasan is used; with libasan6 it has been observed failure for the allocator to return NULL. So we tolerate this."
			exit 0
		else
			exit 1
		fi
	else
		exit 1
	fi
fi
