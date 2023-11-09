AC_PREREQ([2.69])
AC_INIT([librsbpp],[0.0],[michelemartone_AT_users_DOT_sourceforge_DOT_net])
AC_COPYRIGHT([Copyright (c) 2020-2022, Michele Martone])
AM_INIT_AUTOMAKE
LT_INIT
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_CC
AC_PROG_CXX
AC_LANG([C++])
AC_CHECK_HEADERS([vector])
AC_CHECK_HEADERS([iomanip])
if test x"${ac_cv_header_iomanip}" != x"yes" ; then
	AC_MSG_ERROR([It seems like the compiler is missing part of the Standard Library!]);
fi
if test x"${ac_cv_header_vector}" != x"yes" ; then
	AC_MSG_ERROR([It seems like the compiler is missing part of the Standard Template Library!]);
fi
AC_CHECK_HEADERS([gmock/gmock.h])
AC_CHECK_HEADERS([gtest/gtest.h])
AM_EXTRA_RECURSIVE_TARGETS([realclean])
AC_CONFIG_HEADERS([config.h])
AC_C_BIGENDIAN()
AC_OPENMP()
if test x != x"${ac_cv_env_OPENMP_CXXFLAGS_value}" ; then
	OPENMP_CXXFLAGS="${ac_cv_env_OPENMP_CXXFLAGS_value}";
	AC_MSG_NOTICE([Using user-provided OPENMP_CXXFLAGS=${OPENMP_CXXFLAGS}])
fi
AC_ARG_ENABLE(debug, AS_HELP_STRING([--enable-debug],[Compile with debug flags and enable assertions and other internals.]), [if test "x$enableval" = xyes; then enable_debug=yes; else enable_debug=no; fi],[enable_debug=no])
AC_ARG_ENABLE(minimal-librsbpp, AS_HELP_STRING([--enable-minimal-librsbpp],[Enable only minimal librsbpp(EXPERIMENTAL).]), [if test "x$enableval" = xyes; then want_mini_rsbpp="1"; else want_mini_rsbpp="0";fi],[want_mini_rsbpp="0"])
AC_CHECK_PROGS([M4],[gm4 m4 /opt/freeware/bin/m4],[false],[])
AC_ARG_VAR([LIBRSBPP__CPPFLAGS],[Reserved for use by librsb subpackages: do not set this.])
AC_ARG_VAR([LIBRSB_CONFIG],[Path to a librsb-config of a valid librsb installation. Set to e.g. 'false' to force failure, or empty for auto detection.])
if test x"$LIBRSBPP__CPPFLAGS" != x ; then
	export LIBRSB_CONFIG=false;
	AC_MSG_NOTICE([Using internal, configure-set LIBRSBPP__CPPFLAGS variable to ignore LIBRSB_CONFIG.])
	AC_MSG_NOTICE([Also appending LIBRSBPP__CPPFLAGS="${LIBRSBPP__CPPFLAGS}" to CPPFLAGS.]);
	CPPFLAGS="${CPPFLAGS} ${LIBRSBPP__CPPFLAGS}";
fi
if test x"$LIBRSB_CONFIG" = "x" ; then
	AC_PATH_PROG([LIBRSB_CONFIG],[librsb-config],[no])
fi
if test -n "$LIBRSB_CONFIG" -a "x$LIBRSB_CONFIG" != x"no" -a "x$LIBRSB_CONFIG" != x"false" && LIBRSB_CONFIG_VERSION=$("$LIBRSB_CONFIG" --version) ; then
	CPPFLAGS="${CPPFLAGS} ` $LIBRSB_CONFIG --I_opts`"
	CXXFLAGS="${CXXFLAGS} ` $LIBRSB_CONFIG --I_opts`"
dnl	LDFLAGS="${LDFLAGS} `$LIBRSB_CONFIG --ldflags`"
dnl	LIBS="${LIBS} `$LIBRSB_CONFIG --libs --extra_libs`"
	AC_MSG_NOTICE([Found a usable librsb-config ($LIBRSB_CONFIG), using CPPFLAGS=${CPPFLAGS}  CXXFLAGS=${CXXFLAGS}.]);
elif test x"$LIBRSBPP__CPPFLAGS" = x ; then
	AC_MSG_NOTICE([No usable librsb-config found in the path or LIBRSB_CONFIG environment variable.]);
fi
AC_CHECK_HEADERS([rsb.h], [], [], [[
#define RSB_WANT_NO_RSB_TYPES_H 1
]])
if test x"$ac_cv_header_rsb_h" = x"yes" ; then
	CPPFLAGS="${CPPFLAGS} -DRSBPP_HAS_RSB_H=1"
	CXXFLAGS="${CXXFLAGS} -DRSBPP_HAS_RSB_H=1"
else
	if test x"$LIBRSBPP__CPPFLAGS" = x ; then
		AC_MSG_WARN([No usable rsb.h found: using own (default) defines.]);
	else
		AC_MSG_ERROR([No usable rsb.h found, despite LIBRSBPP__CPPFLAGS internal flag ?!]);
	fi
fi
AC_CHECK_HEADERS([getopt.h])
AC_CHECK_HEADERS([unistd.h])
AC_CHECK_HEADERS([omp.h])
AC_CHECK_FUNC([rsb_lib_init], AC_DEFINE([HAVE_LIBRSB], [1], [have librsb]))
AC_CHECK_LIB([pthread],[pthread_cancel],[LIBS="${LIBS} -lpthread";])
AC_CHECK_LIB([gtest],[main],[true]) dnl with g++-10 one may look for _ZN7testing4TestD0Ev but this is mangled and therefore volatile.
AC_CHECK_LIB([gtest_main],[main],[true])
AC_CHECK_LIB([stdc++fs],[main],[LIBS="$LIBS -lstdc++fs";]) dnl With g++-10 one may look for _ZNSt10filesystem9file_sizeERKNS_7__cxx114pathE but this is mangled and therefore volatile.
AC_CHECK_HEADERS([filesystem], AC_DEFINE([HAVE_FILESYSTEM], [1], [Define to 1 if have <filesystem>]))
dnl AC_CHECK_LIB([stdc++],[_ZNKSs13get_allocatorEv],[LIBS="${LIBS} -lstdc++";])
if test x"$ac_cv_lib_gtest_main" = x"yes" -a x"$ac_cv_lib_gtest_main_main" = x"yes" -a x"$ac_cv_header_gmock_gmock_h" = x"yes" -a x"$ac_cv_header_gtest_gtest_h" = x"yes"; then
	rsb_have_gtest=1;
	RPPGT_LDADD="-lgtest_main -lgtest -lpthread"
	AC_MSG_NOTICE([Will use Google Test.])
else
	AC_MSG_NOTICE([Will not use Google Test.])
	rsb_have_gtest=0;
fi
AC_CACHE_CHECK([whether you have std::thread],
        [have_cv_std_thread], AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <thread>
int main(const int argc, char ** argv) { std::thread* tvp; return 0; }
]])],[have_cv_std_thread=yes],[have_cv_std_thread=no],[have_cv_std_thread=no]))
if test X"$have_cv_std_thread" = X"yes"; then
	AC_DEFINE([RSBPP_HAVE_STD_THREAD],[1],[Have std::thread.])
else
	AC_MSG_NOTICE([No std::thread found.])
fi
AC_CACHE_CHECK([whether you have std::mutex],
        [have_cv_std_mutex], AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <mutex>
int main(const int argc, char ** argv) { std::mutex* tvp; return 0; }
]])],[have_cv_std_mutex=yes],[have_cv_std_mutex=no],[have_cv_std_mutex=no]))
if test X"$have_cv_std_mutex" = X"yes"; then
	AC_DEFINE([RSBPP_HAVE_STD_MUTEX],[1],[Have std::mutex.])
else
	AC_MSG_NOTICE([No std::mutex found.])
fi
AH_TEMPLATE([RSBPP_PRESET_SPMM_BLOCKINGS])
AC_ARG_ENABLE(spmm-unrolls, AS_HELP_STRING([--enable-spmm-unrolls],[Unroll SpMM kernels (default).]), [want_spmm_unrolls="$enableval";], [want_spmm_unrolls="yes";])
AH_TEMPLATE([RSB_HAVE_GTEST])
AC_ARG_ENABLE(programs, AS_HELP_STRING([--enable-programs],[Build test and example programs.]), [want_rsb_programs="$enableval";], [want_rsb_programs="yes";])
if test "x$want_rsb_programs" = xno; then
	AC_MSG_NOTICE(Will not build any programs.)
fi
if test x$want_spmm_unrolls = x"yes"; then
	rsbpp_preset_spmm_blockings=4;
else
	rsbpp_preset_spmm_blockings=1;
fi
if test "x$want_mini_rsbpp" = x"1"; then
	AC_MSG_NOTICE([Enabling librsbpp library only (no extras).])
	CXXFLAGS="${CXXFLAGS} -DRSBPP_WANT_MINIMAL_LIBRSBPP=1"
fi
AC_SUBST(RPPGT_LDADD,"$RPPGT_LDADD")
dnl AC_DEFINE_UNQUOTED([RSBPP_PRESET_SPMM_BLOCKINGS],[${rsbpp_preset_spmm_blockings}],[])
AC_SUBST([RSBPP_PRESET_SPMM_BLOCKINGS],["$rsbpp_preset_spmm_blockings"],[])
AM_CONDITIONAL([HAVE_WHICH],[which ${SHELL}])
AM_CONDITIONAL([WANT_NO_BLOCKINGS],[test x"${rsbpp_preset_spmm_blockings}" == x"1"])
AC_DEFINE_UNQUOTED([RSB_HAVE_GTEST],[${rsb_have_gtest}],[])
AM_CONDITIONAL([HAVE_GTEST],[test x"${rsb_have_gtest}" = x"1"])
AM_CONDITIONAL([WANT_MINIMAL_LIBRSBPP],[test x"${want_mini_rsbpp}" == x"1"])
AM_CONDITIONAL([WANT_PROGRAMS],[test x"${want_rsb_programs}" != x"no"])
AM_CONDITIONAL([WANT_DBG],[test x"${enable_debug}" == x"yes" -a x"$ac_cv_cxx_compiler_gnu" = x"yes"])
AM_CONDITIONAL([HAVE_LIBRSB_M4_MACROS_DIR],[test x"$M4" != xfalse -a -f ../wisdom.m4 -a -f ../rsb_misc.m4 -a -f ../do_unroll.m4])
AM_CONDITIONAL([WANT_CONVENIENCE_ONLY],[test x"$LIBRSBPP__CPPFLAGS" != x ])
if test x"$host_cpu" == x"$build_cpu"; then
	ac_cv_cross_compiling="no";
else
	ac_cv_cross_compiling="yes";
fi
AM_CONDITIONAL([HAVE_RUNNABLE_TESTS],[test x"$ac_cv_header_getopt_h" = x"yes" -a x"$ac_cv_cross_compiling" != x"yes" ])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
AC_MSG_NOTICE([Successfully created a Makefile.])
if test x"$M4" = x"false"; then
	AC_MSG_WARN([You do not have an m4 processor. As a consequence, code-generation is disabled.])
fi
