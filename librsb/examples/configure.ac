dnl Example minimal configure.ac you can use with your LIBRSB-based program.
AC_PREREQ([2.69])
AC_INIT([examples],[0.0],[michele martone])
AC_COPYRIGHT([Copyright (c) 2021-2022, Michele Martone])
AM_INIT_AUTOMAKE()
AC_OPENMP()
AC_PROG_CC()
AC_PROG_FC()
AC_LANG(C++)
AC_PROG_CXX()
AC_CANONICAL_BUILD()
AC_CANONICAL_HOST()
AC_CONFIG_HEADERS([config.h])
AC_ARG_VAR([LIBRSB_CONFIG],[Path to a librsb-config of a valid librsb installation.])
AC_CHECK_PROG([HAVE_LIBRSB_CONFIG],[librsb-config],[yes],[no])
AC_PATH_PROG([LIBRSB_CONFIG],[librsb-config],[no])
if test x$LIBRSB_CONFIG == x"no" ; then
	AC_MSG_ERROR([Where is librsb-config ?]);
else
	AC_MSG_NOTICE([Adding LIBRSB_CONFIG --I_opts/--ldflags/--libs --extra_libs  to CPPFLAGS/LDFLAGS/LIBS.])
	CPPFLAGS="${CPPFLAGS} `${LIBRSB_CONFIG} --I_opts`"
	CXXFLAGS="${CXXFLAGS} `${LIBRSB_CONFIG} --cxxflags`"
	LDFLAGS="${LDFLAGS} `${LIBRSB_CONFIG} --ldflags`"
	LIBS="${LIBS} `${LIBRSB_CONFIG} --libs --extra_libs`"
fi
AC_CHECK_HEADERS([rsb.h], [break], AC_MSG_ERROR([No rsb.h ?!]))
AC_CHECK_HEADERS([rsb.hpp], [break], AC_MSG_NOTICE([No rsb.hpp.]))
AC_CHECK_FUNC([rsb_lib_init], AC_DEFINE([HAVE_LIBRSB], [1], [have librsb]), AC_MSG_ERROR([No usable librsb ?!]))
AM_CONDITIONAL([WANT_CROSS_COMPILE],[test x"$build_cpu" != x"host_cpu" ])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
AC_MSG_NOTICE([Created a Makefile.])
