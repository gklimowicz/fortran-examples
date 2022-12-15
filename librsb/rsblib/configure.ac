AC_PREREQ([2.69])
AC_INIT([rsblib],[0.0],[michelemartone_AT_users_DOT_sourceforge_DOT_net])
AC_COPYRIGHT([Copyright (c) 2017-2022, Michele Martone])
AM_INIT_AUTOMAKE
LT_INIT
AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
AC_PROG_CXX
AC_OPENMP()
if test x != x"${ac_cv_env_OPENMP_CXXFLAGS_value}" ; then
	OPENMP_CXXFLAGS="${ac_cv_env_OPENMP_CXXFLAGS_value}";
	AC_MSG_NOTICE([Using user-provided OPENMP_CXXFLAGS=$OPENMP_CXXFLAGS])
fi
AC_CHECK_HEADERS([gmock/gmock.h])
AC_CHECK_HEADERS([gtest/gtest.h])
AC_CONFIG_HEADERS([config.h])
AC_C_BIGENDIAN()
AC_PATH_PROG([DOXYGEN],[doxygen],[no])
AC_CHECK_PROGS([MAKE],[gmake make /opt/freeware/bin/make],[false],[])
AC_ARG_VAR([LIBRSB__BASE_DIR],[Reserved for use by librsb subpackages: do not set this.])
if test -n "$LIBRSB__BASE_DIR"; then
	AC_MSG_NOTICE([Assuming building as subdir of librsb, and that librsb will be built later.])
	want_rsblib="yes";
fi
BPWD="`cd ..;pwd`";
AC_ARG_WITH(rsblib, AS_HELP_STRING([--with-rsblib],[Link to librsb build dir when building in a subdir (EXPERIMENTAL, developer only).]), [if test "x$withval" != xno; then want_rsblib="$withval" ; else want_rsblib="no"; fi ], [if test x"$want_rsblib" = x; then want_rsblib="no"; fi;])
if test x"$want_rsblib" != x"no" ; then
	LIBRSB_LOCAL_LIBS="${BPWD}/librsb.la";
	CPPFLAGS="${CPPFLAGS} -I${BPWD}"
	ac_cv_func_rsb_lib_init="yes";
else
	dnl
	AC_ARG_VAR([LIBRSB_CONFIG],[Path to a librsb-config of a valid librsb installation.])
	AC_CHECK_PROG([HAVE_LIBRSB_CONFIG],[librsb-config],[yes],[no])
	AC_PATH_PROG([LIBRSB_CONFIG],[librsb-config],[no])
	if test x$LIBRSB_CONFIG == x"no" ; then
		AC_MSG_NOTICE([Where is librsb-config ?]);
	else
		want_rsblib=yes;
		CPPFLAGS="${CPPFLAGS} `$LIBRSB_CONFIG --I_opts`"
		LDFLAGS="${LDFLAGS} `$LIBRSB_CONFIG --ldflags`"
		LIBS="${LIBS} `$LIBRSB_CONFIG --libs --extra_libs`"
	fi
	AC_CHECK_HEADERS([rsb.h], [], [], [[
	#define RSB_WANT_NO_RSB_TYPES_H 1
	]])
	AC_CHECK_FUNC([rsb_lib_init], AC_DEFINE([HAVE_LIBRSB], [1], [have librsb]))
fi
AC_CHECK_HEADERS([getopt.h])
AC_CHECK_HEADERS([unistd.h])
AC_CHECK_HEADERS([omp.h])
AC_CHECK_LIB([pthread],[pthread_cancel],[LIBS="${LIBS} -lpthread";])
AC_CHECK_LIB([gtest],[main],[true]) dnl with g++-10 one may look for _ZN7testing4TestD0Ev but this is mangled and therefore volatile.
AC_CHECK_LIB([gtest_main],[main],[true])
AC_CHECK_HEADERS([filesystem], AC_DEFINE([HAVE_FILESYSTEM], [1], [Define to 1 if have <filesystem>]))
AH_TEMPLATE([RSB_HAVE_GTEST])
AC_ARG_ENABLE(debug, AS_HELP_STRING([--enable-debug],[Compile with debug flags and enable assertions and other internals.]), [if test "x$enableval" = xyes; then enable_debug=yes; else enable_debug=no; fi],[enable_debug=no])
AC_ARG_ENABLE(matrix-types, AS_HELP_STRING([--enable-matrix-types],[Alias of the librsb equivalent option.]), [want_matrix_types="$enableval"],[want_matrix_types=""])
AC_ARG_ENABLE(programs, AS_HELP_STRING([--enable-programs],[Build test and example programs.]), [want_rsb_programs="$enableval";], [want_rsb_programs="yes";])
if test x"$ac_cv_func_rsb_lib_init" = x"no" ; then
	AC_MSG_NOTICE([No linkable librsb library detected.]);
	want_rsb_programs=no;
fi
if test "x$want_rsb_programs" = xno; then
	AC_MSG_NOTICE(Will not build any programs.)
fi
dnl
AC_CHECK_LIB([stdc++fs],[main],[true;]) dnl With g++-10 one may look for _ZNSt10filesystem9file_sizeERKNS_7__cxx114pathE but this is mangled and therefore volatile.
if test x"$ac_cv_lib_gtest_main" = x"yes" -a x"$ac_cv_lib_gtest_main_main" = x"yes" -a x"$ac_cv_header_gmock_gmock_h" = x"yes" -a x"$ac_cv_header_gtest_gtest_h" = x"yes"; then
	rsb_have_gtest=1;
	RGT_LDADD="-lgtest_main -lgtest -lpthread"
	if test x"$ac_cv_lib_stdcppfs_main" = x"yes" ; then
		RGT_LDADD="$RGT_LDADD -lstdc++fs";
	fi
else
	rsb_have_gtest=0;
fi
if test x"$want_rsblib" != x"no" ; then
	# adding this now: assuming configure called *before* librsb.la built.
	LIBS="${LIBS} ${LIBRSB_LOCAL_LIBS}"
	AC_MSG_NOTICE([Assuming you are yet to build librsb.la. (set LIBS=$LIBS)])
fi
if test x"$ac_cv_header_getopt_h" != x"yes" ; then
	AC_MSG_WARN([No <getopt.h> header detected. Do NOT expect the test suite to work!])
fi
AC_ARG_ENABLE(c++-examples, AS_HELP_STRING([--disable-c++-examples],[C++ example programs building.]), [if test "x$enableval" = xno; then enable_cpp_examples=no; else enable_cpp_examples=yes; fi;],[enable_cpp_examples=yes;])
AC_SUBST(RGT_LDADD,"$RGT_LDADD")
AC_DEFINE_UNQUOTED([RSB_HAVE_GTEST],[${rsb_have_gtest}],[])
AM_CONDITIONAL([HAVE_WHICH],[which ${SHELL}])
if test x"$host_cpu" == x"$build_cpu"; then
	ac_cv_cross_compiling="no";
else
	ac_cv_cross_compiling="yes";
fi
AM_CONDITIONAL([HAVE_RUNNABLE_TESTS],[test x"$ac_cv_header_getopt_h" = x"yes" -a x"$ac_cv_cross_compiling" != x"yes" ])
AM_CONDITIONAL([HAVE_STDCPPFS],[test x"${ac_cv_lib_stdcppfs_main}" = x"yes"])
AM_CONDITIONAL([HAVE_GTEST],[test x"${rsb_have_gtest}" = x"1"])
AM_CONDITIONAL([WANT_PROGRAMS],[test x"${want_rsb_programs}" != x"no"])
AM_CONDITIONAL([WANT_EXAMPLES],[test x"${enable_cpp_examples}" != x"no" -a x"${want_rsb_programs}" != x"no"])
AM_CONDITIONAL([WANT_CUSTOM_NUMERICAL_TYPES],[test x"${want_matrix_types}" != x""])
AM_CONDITIONAL([HAVE_DOXYGEN],[test x"${DOXYGEN}" != x"no"])
AM_CONDITIONAL([WANT_RSBLIB_AS_SUBDIR],[test x"$want_rsblib" != x"no"])
AM_CONDITIONAL([WANT_DBG],[test x"${enable_debug}" == x"yes" -a x"$ac_cv_cxx_compiler_gnu" = x"yes"])
AM_EXTRA_RECURSIVE_TARGETS([realclean])
AC_PROG_MKDIR_P()
AC_CONFIG_FILES([Makefile examples/Makefile])
AC_OUTPUT
AC_MSG_NOTICE([Created a Makefile.])
