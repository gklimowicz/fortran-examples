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

# Portability notes:
#  avoid E\+, use EE*
#  avoid \n, use a newline
#  avoid \s, use [[:space:]]
#  avoid [A-F0-9], use [[:xdigit:]]
#  avoid \<E\>, use context

set -e
set -o pipefail
SRCDIR=
BLDDIR=
if test $# = 0 ; then SRCDIR=. ; else SRCDIR="$1"; fi
if test $# = 1 ; then BLDDIR=. ; else BLDDIR="$2"; fi
IF=${SRCDIR}/rsb.h
TF=${BLDDIR}/rsb_types.h
CF=${BLDDIR}/rsb-config.h
CH2ICFB=./ch2icfb
test -x ${CH2ICFB} || { echo " [!] ${0}: missing executable ${CH2ICFB}!" 1>&2 ; false; }
test -f ${TF} || { echo " [!] ${0}: missing generated header ${TF}!" 1>&2 ; false; }
SED=${SED:=sed}
GREP=${GREP:=grep}
(
${CH2ICFB} < ${IF} | grep -v 'END MODULE rsb' 
grep rsb_blas_file_mtx_load ${SRCDIR}/rsb_libspblas.h | CH2ICFB_NH=1 ${CH2ICFB}
SHEXP='s/0x\([[:xdigit:]][[:xdigit:]]*\)/INT(Z"0\1",C_INT)/g'
IPD='INTEGER(C_INT),PARAMETER::'
IPD2='INTEGER(C_INT),PARAMETER::'
FD='s/^\([[:graph:]][[:graph:]]*\) \([[:graph:]][[:graph:]]*\)/'"${IPD2}"'\1=\2/g'
CLEANUP='s/[[:space:]][[:space:]]*/ /g;s/^,//g;s/=//g;s/\/.*$//g;s/^[[:space:]]*//g;s/#define *//g'
D2N='s/#define //g'
DS='^#define '
SEE='s/\(PARAMETER::*\) *\(RSB[A-Z_0-9]*\)\(.*$\)/\1\2\3 !< See #\2./g'
IC='      '
SHORTEN_DC='s/\(::\)/\&\
'"${IC}"'\&\1/g;'
SHORTEN_EX='s/\([A-Z_]++*\)/\1\&\
'"${IC}"'\&/g;'
SHORTEN_PA='s/\( *:: *[A-Z_][A-Z_]*\)/\1\&\
'"${IC}"'\&/g;'"$SHORTEN_EX""${SHORTEN_DC}"
#SHORTEN_PM='s/\([=+]\)/\&\n'"${IC}"'\&\1/g;'
SHORTEN_TK='s/[[:space:]][[:space:]]*/\&\
'"${IC}"'\&/g;'
NOTS='s/[[:space:]]*$//g;'

echo '! Error values '
${SED} 's/[[:space:]][[:space:]]*/ /g;s/^\(.define\) \(RSB_ERR[^ ]*\) RSB_ERR_CAST(0x\([^ ]*\))$/DEFINE \2 = -INT(Z"0\3",C_INT)/g;s/RSB_ERR_CAST/-/g;s/DEFINE */'"${IPD}"'/g;' < ${IF} | ${GREP} '^ *INTE.*RSB_ERR' | ${GREP} -v 'RSB_ERRS_UNSUPPORTED_FEATURES\|RSB_ERR_TO_PROGRAM_ERROR'  | ${SED} "${SEE}"| ${SED} "${NOTS}${SHORTEN_PA}"

echo '! Matrix flags values '
${GREP} RSB_FLAG_ ${IF} | ${GREP} -v '\\$' | ${GREP} '^.define' | ${SED} 's/[[:space:]][[:space:]]*/ /g;'"${SHEXP}" | ${GREP} -v '\/.*' | ${SED} 's/[[:space:]][[:space:]]*/ /g;s/^\(.define\) \(RSB_FLAG[^[:space:]]*\) \(INT(Z[^[:space:]]*\)$/DEFINE\2 = \3/g;s/DEFINE/'"${IPD}"'/g;'  | ${GREP} '^ *INTE.*RSB_FLAG'  | ${SED} "${SEE}"| ${SED} "${NOTS}${SHORTEN_PA}"

echo '! Composite flags '
${GREP} RSB_FLAG_ ${IF} | ${GREP} -v '[ 	]0x'  | ${SED} 's/[[:space:]][[:space:]]*/ /g;s/|/+/g;s/^\(.define\)[[:space:]]\(RSB_FLAG[^	 ]*\)[[:space:]]\(.*$\)/DEFINE \2 = \3/g;s/^ *//g;s/DEFINE/'"${IPD2}"'/g' | ${GREP} '^ *INTE.*RSB_FLAG' | ${SED} "${SEE}" | ${SED} "${NOTS}${SHORTEN_PA}"

echo '! Transposition constants '
${GREP} "${DS} *"'RSB_TRANSPOSITION_[NTC]' "${TF}" | ${SED} "${CLEANUP};${D2N};${SHEXP};${FD}"

echo '! Numerical types constants '
${GREP} "${DS} *"'RSB_NUMERICAL_TYPE_FORTRAN_' "${TF}" | ${SED} "${CLEANUP};${D2N};${SHEXP};${FD};s/_FORTRAN//g" | ${SED} 's/C_INT/C_SIGNED_CHAR/g' | ${SED} "${SHORTEN_DC}"
# ( ${GREP} "${DS} *"'RSB_WANT_LONG_IDX 1' "${CF}" || echo '#define RSB_WANT_LONG_IDX 0' ) | ${SED} 's/RSB_WANT_LONG_IDX/RSB_IDX_KIND/g;s/1/8/g;s/0/4/g' | ${SED} "${CLEANUP};${D2N};${SHEXP};${FD};s/_FORTRAN//g" | ${SED} 's/C_INT/C_SIGNED_CHAR/g' | ${SED} "${SHORTEN_DC}"

echo '! Other enumerations constants '
${GREP} -E '^(.define|[ ,]*) *RSB_(IO_WANT|MARF|PRECF|EXTF|MIF|ELOPF)_' ${IF} | ${SED} "${CLEANUP};${SHEXP};${FD};${SEE}"| ${SED} "${NOTS}${SHORTEN_PA}${SHORTEN_EX}"
${GREP} -E '^(.define) *RSB_(NULL)_' ${IF} | ${SED} "${CLEANUP};${SHEXP};${FD};${SEE}" | ${SED} "${NOTS}${SHORTEN_PA}"| ${SED} 's/=NULL/=C_NULL_PTR/g;s/INTEGER(C_INT)/TYPE(C_PTR)/g'

echo 'END MODULE rsb'
) | ${SED} 's/^/      /g;s/^\( *!\)/!/g;s/^[[:space:]]*#/#/g'
