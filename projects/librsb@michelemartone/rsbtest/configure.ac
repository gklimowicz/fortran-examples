AC_PREREQ([2.69])
AC_INIT([rsbtest],[0.0],[michelemartone_AT_users_DOT_sourceforge_DOT_net])
AC_COPYRIGHT([Copyright (c) 2017-2022, Michele Martone])
AM_INIT_AUTOMAKE
LT_INIT
AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
AC_PROG_CXX
AC_OPENMP
if test x != x"${ac_cv_env_OPENMP_CXXFLAGS_value}" ; then
	OPENMP_CXXFLAGS="${ac_cv_env_OPENMP_CXXFLAGS_value}";
	AC_MSG_NOTICE([Using user-provided OPENMP_CXXFLAGS=$OPENMP_CXXFLAGS])
fi
AC_CONFIG_HEADERS([config.h])
AC_C_BIGENDIAN()
dnl AC_CHECK_TOOL(PDFLATEX, [pdflatex], [false])
AC_CHECK_TOOL(LATEX, [latex], [false])
AC_CACHE_CHECK([whether you have C++17],
        [have_cv_cxx17], AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#if (__cplusplus >= 201703L)
int main() { }
#else
#error No C++17 supported, it seems!
#endif
]])],[have_cv_cxx17=yes],[have_cv_cxx17=no],[have_cv_cxx17=no]))
if test X$have_cv_cxx17 = Xyes; then
	true
else
	AC_MSG_WARN([No C++17 support detected.])
fi
AC_ARG_VAR([LIBRSB__BASE_DIR],[Reserved for use by librsb subpackages: do not set this.])
BPWD="`cd ..;pwd`";
CPPFLAGS="${RSBTEST_CPPFLAGS} ${CPPFLAGS}"
AC_ARG_WITH(librsbpp, AS_HELP_STRING([--with-librsbpp],[Link to librsb/librsbpp build dir when building in a subdir (EXPERIMENTAL, developer only).]), [if test "x$withval" != xno; then want_librsbpp="$withval" ; else want_librsbpp="no"; fi ], [want_librsbpp="no"])
if test x"$want_librsbpp" != x"no" ; then
	CPPFLAGS="${CPPFLAGS} -I${BPWD}/librsbpp/ "
	AC_MSG_NOTICE([Updated CPPFLAGS to have librsbpp dir: ${CPPFLAGS}.])
fi
if test -n "$LIBRSB__BASE_DIR"; then
	want_rsbtest="yes";
fi
AC_ARG_WITH(rsbtest, AS_HELP_STRING([--with-rsbtest],[Link to librsb build dir when building in a subdir (EXPERIMENTAL, developer only).]), [if test "x$withval" != xno; then want_rsbtest="$withval" ; else want_rsbtest="no"; fi ], [if test "x$want_rsbtest" = x; then want_rsbtest="no"; fi;])
if test x"$want_rsbtest" != x"no" ; then
	AC_SUBST([LIBRSB_LOCAL_LIBS],"${BPWD}/librsb.la")
	AC_MSG_NOTICE([Assuming building as subdir of librsb: link to $LIBRSB_LOCAL_LIBS.])
	AM_CONDITIONAL([WITH_RSBTEST],[true])
else
	AM_CONDITIONAL([WITH_RSBTEST],[false])
	AC_SUBST([LIBRSB_LOCAL_LIBS],"")
	dnl
	AC_ARG_VAR([LIBRSB_CONFIG],[Path to a librsb-config of a valid librsb installation.])
	AC_ARG_VAR([RSBTEST_CPPFLAGS],[Specific CPPFLAGS.])
	AC_CHECK_PROG([HAVE_LIBRSB_CONFIG],[librsb-config],[yes],[no])
	AC_PATH_PROG([LIBRSB_CONFIG],[librsb-config],[no])
	if test x$LIBRSB_CONFIG == x"no" ; then
		AC_MSG_ERROR([Where is librsb-config ?]);
	else
		AC_MSG_NOTICE([Adding LIBRSB_CONFIG --I_opts/--ldflags/--libs --extra_libs  to CPPFLAGS/LDFLAGS/LIBS.])
		CPPFLAGS="${CPPFLAGS} `${LIBRSB_CONFIG} --I_opts`"
		LDFLAGS="${LDFLAGS} `${LIBRSB_CONFIG} --ldflags`"
		LIBS="${LIBS} `${LIBRSB_CONFIG} --libs --extra_libs`"
	fi
fi
AC_ARG_ENABLE(long-indices, AS_HELP_STRING([--enable-long-indices],[Enable long (64 bit) indices for rsb_coo_idx_t and rsb_nnz_idx_t (EXPERIMENTAL).]), [want_long_idx="1"],[want_long_idx="0"])
if test "x$want_long_idx" = x"1"; then
	AC_MSG_NOTICE([Enabling long types flags for rsb_coo_idx_t and rsb_nnz_idx_t.])
fi
AC_ARG_ENABLE(preset, AS_HELP_STRING([--enable-preset],[Compile with preset options (developer only)]), [ if test x$enableval = xyes; then 
	AC_MSG_NOTICE([Will compile with preset benchmarking options.]);
		CPPFLAGS="${CPPFLAGS} -std=c++17 -O3 -fopt-info-vec=vec.opt  -Wno-sign-compare" ;
		if test -n "${MKL_INC}"; then CPPFLAGS="${CPPFLAGS} ${MKL_INC}"; fi
		RSBLIBDIR="${HOME}/src/rsblib/"
		if test -d "${RSBLIBDIR}" ; then CPPFLAGS="${CPPFLAGS} -I${RSBLIBDIR}"; fi
		LIBRSBPPDIR="${HOME}/src/librsbpp/"
		if test -d "${LIBRSBPPDIR}" ; then CPPFLAGS="${CPPFLAGS} -I${LIBRSBPPDIR}"; fi
		#if test -n "${MKL_LIB}"; then LIBS="${LIBS} ${MKL_LIB}"; fi
		if test -n "${MKL_LIBDIR}"; then LIBS="${LIBS} -Wl,--start-group ${MKL_LIBDIR}/libmkl_intel_lp64.a ${MKL_LIBDIR}/libmkl_core.a ${MKL_LIBDIR}/libmkl_intel_thread.a -Wl,--end-group -L ${MKL_LIBDIR}/../../../compiler/lib/intel64 -liomp5 -lpthread -lm"; fi
fi ])
AC_ARG_WITH(mkl, AS_HELP_STRING([--with-mkl],[Specify if MKL shall be used if detected (--without-mkl to disable). ]), [if test "x$withval" = xno; then want_mkl_libs= ; else if test "x$withval" = xyes; then want_mkl_libs="yes" ; else want_mkl_libs="$withval" ; fi  ; fi], [true] )
AC_CHECK_LIB([stdc++fs],[main],[LIBS="${LIBS} -lstdc++fs"])
dnl AC_CHECK_LIB([stdc++fs],[_ZNSt12experimental10filesystem2v17__cxx1116filesystem_errorD1Ev],[LIBS="${LIBS} -lstdc++fs"])
	AC_CHECK_HEADERS([mkl/mkl.h])
	AC_CHECK_HEADERS([mkl/mkl_spblas.h])
if test x"$want_mkl_libs" = x"yes" ; then
	dnl AC_CHECK_LIB([stdc++fs],[std::experimental::filesystem::is_directory])
	RSB_MKL_LIBS_GNU64_SET=" -fopenmp -lpthread -Wl,--start-group,-lmkl_intel_lp64,-lmkl_gnu_thread,-lmkl_core,--end-group -lm"
	AH_TEMPLATE([HAVE_MKL_GET_VERSION])
	AC_CHECK_LIB([mkl_core], [mkl_get_version], [LIBS="${LIBS} ${RSB_MKL_LIBS_GNU64_SET}";] AC_DEFINE([HAVE_MKL_GET_VERSION], [1], []), [LIBS="${LIBS}";] , ["${RSB_MKL_LIBS_GNU64_SET}"])
	AC_DEFINE([WANT_MKL], [1], [want MKL])
fi
AM_CONDITIONAL([USE_MKL_SPBLAS],[test x"${want_mkl_lib}" == x"yes" -a x"$ac_cv_header_mkl_mkl_spblas_h" == x"yes"])
if test x"$want_mkl_libs" = x"yes" -a x"$want_long_idx" = x"1"; then
	AC_MSG_ERROR([Sorry, as of now, you can't use MKL and long librsb indices!]);
fi
AC_ARG_ENABLE(minimal-librsbpp, AS_HELP_STRING([--enable-minimal-librsbpp],[Enable only minimal librsbpp (EXPERIMENTAL).]), [if test "x$enableval" = xyes; then want_mini_rsbpp="1"; else want_mini_rsbpp="0";fi],[want_mini_rsbpp="0"])
if test "x$want_mini_rsbpp" = x"1"; then
	AC_MSG_NOTICE([Enabling librsbpp library only usage (no extras).])
	CXXFLAGS="${CXXFLAGS} -DRSBT_WANT_MINIMAL_LIBRSBPP=1"
fi
AM_CONDITIONAL([WANT_MINIMAL_LIBRSBPP],[test x"${want_mini_rsbpp}" == x"1"])
AC_CHECK_HEADERS([rsb.h], [], [], [[
#define RSB_WANT_NO_RSB_TYPES_H 1
]])
AC_CHECK_HEADERS([rsb.hpp], [], [], [[
#define RSB_RSB_HPP_INCLUDED 1	/* This checks for rsb.hpp presence only -- included rsb_types.h may as well not exist yet at configure time. */
]])
AC_ARG_ENABLE(matrix-types, AS_HELP_STRING([--enable-matrix-types],[Alias of the librsb equivalent option.]), [want_matrix_types="$enableval"],[want_matrix_types=""])
AC_ARG_ENABLE(programs, AS_HELP_STRING([--enable-programs],[Install rsbtest.]), [want_rsb_programs="$enableval";], [want_rsb_programs="no";])
AM_CONDITIONAL([WANT_PROGRAMS],[test x"${want_rsb_programs}" != x"no"])
AM_CONDITIONAL([WANT_CUSTOM_NUMERICAL_TYPES],[test x"${want_matrix_types}" != x""])
AC_CHECK_HEADERS([rsbpp.hpp])
AC_CHECK_HEADERS([getopt.h])
AC_CHECK_HEADERS([sys/utsname.h])
AC_CHECK_HEADERS([omp.h])
AC_CHECK_HEADERS([mutex])
AC_CHECK_HEADERS([filesystem], AC_DEFINE([HAVE_FILESYSTEM], [1], [Define to 1 if have <filesystem>]))
AC_CHECK_FUNC([rsb_lib_init], AC_DEFINE([HAVE_LIBRSB], [1], [have librsb]))
AM_CONDITIONAL([WITH_CXX17_FILESYSTEM],[test x"$have_cv_cxx17" = x"yes" -a x"$ac_cv_header_filesystem" = x"yes"])
AM_CONDITIONAL([HAVE_RSBPP_HPP],[test x"$ac_cv_header_rsbpp_hpp" = x"yes" -a x"$want_mini_rsbpp" != x"1"])
if test x"$have_cv_cxx17" != x"yes" -a x"$ac_cv_header_rsb_hpp" = x"yes" ; then
	AC_MSG_WARN([You need C++17 or later support to work with <rsb.hpp>! Turning <rsb.hpp> support off.]);
fi
AM_CONDITIONAL([HAVE_RSB_HPP],[test x"$have_cv_cxx17" = x"yes" -a x"$ac_cv_header_rsb_hpp" = x"yes"])
AM_CONDITIONAL([HAVE_DEV_NULL],[test -w /dev/null])
if test x"$ac_cv_header_getopt_h" != x"yes" ; then
	AC_MSG_WARN([No <getopt.h> header detected. Do NOT expect the test suite to work!])
fi
if test x"$host_cpu" == x"$build_cpu"; then
	ac_cv_cross_compiling="no";
else
	ac_cv_cross_compiling="yes";
fi
AM_CONDITIONAL([HAVE_RUNNABLE_TESTS],[test x"$ac_cv_header_getopt_h" = x"yes" -a x"$ac_cv_cross_compiling" != x"yes"])
AC_CHECK_PROG([have_dd],[dd],[yes],[no])
if test x"$ac_cv_prog_have_dd" != x"no" ; then
	if TEST=`echo '' | dd bs=16M status=none iflag=fullblock` ; then
		AM_CONDITIONAL([HAVE_DD],[true])
		AM_CONDITIONAL([HAVE_FULL_DD],[true])
	elif TEST=`echo '' | dd ` ; then
		AM_CONDITIONAL([HAVE_DD],[true])
		AM_CONDITIONAL([HAVE_FULL_DD],[false])
		AC_MSG_WARN([Your dd seems verbose and so will be tests..])
	else
		AM_CONDITIONAL([HAVE_DD],[false])
		AM_CONDITIONAL([HAVE_FULL_DD],[false])
		AC_MSG_WARN([No usable dd program detected!])
	fi
else
	AM_CONDITIONAL([HAVE_DD],[false])
	AM_CONDITIONAL([HAVE_FULL_DD],[false])
	AC_MSG_WARN([No usable dd program detected!])
fi
AM_CONDITIONAL([HAVE_WHICH],[which ${SHELL}])
AM_CONDITIONAL([HAVE_TIMEOUT],[which timeout])
AM_EXTRA_RECURSIVE_TARGETS([realclean])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
AC_MSG_NOTICE([Created a Makefile.])
