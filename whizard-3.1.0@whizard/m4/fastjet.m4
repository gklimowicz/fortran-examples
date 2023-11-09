dnl fastjet.m4 -- checks for fastjet library
dnl Includes the fastjet.m4 code from the fastjet webpage
dnl

AC_DEFUN([WO_PROG_FASTJET],
[dnl
AC_REQUIRE([AC_PROG_CXX])

AC_ARG_ENABLE([fastjet],
  [AS_HELP_STRING([--enable-fastjet],
    [enable FastJet for handling event data [[no]]])],
  [], [enable_fastjet="no"])

if test "$enable_fastjet" = "yes"; then
  ACX_CHECK_FASTJET()
  if test "${fjconfig}" = "no"; then
     enable_fastjet="no"
  else
     wo_require_stdcpp="yes"
     AC_MSG_CHECKING([the FastJet version])
     FASTJET_VERSION=`${fjconfig} --version`
     AC_MSG_RESULT([$FASTJET_VERSION])
     AC_SUBST([FASTJET_VERSION])
  fi
else
  AC_MSG_CHECKING([for FASTJET])
  AC_MSG_RESULT([(disabled)])
fi

FASTJET_AVAILABLE_FLAG=$enable_fastjet

AC_SUBST([FASTJET_CXXFLAGS])
AC_SUBST([FASTJET_LIBS])
AC_SUBST([FASTJET_AVAILABLE_FLAG])

AM_CONDITIONAL([FASTJET_AVAILABLE], [test "$enable_fastjet" = "yes"])
])

dnl #####################################################################
dnl CHECK FASTJET BEGIN
dnl
dnl This script can be used in configure scripts to check for the
dnl usability of the FastJet library.
dnl
dnl By defaults, it searches the FastJet library in standard system
dnl locations but an alternative path can be specified using the
dnl --with-fastjet=... configure option
dnl
dnl If FastJet is found and functional, the variables FASTJET_CXXFLAGS
dnl and FASTJET_LIBS are set
AC_DEFUN([ACX_CHECK_FASTJET],
[
dnl ckeck if a directory is specified for FastJet
AC_ARG_WITH(fastjet,
            [AS_HELP_STRING([--with-fastjet=dir],
                            [assume the given directory for FastJet])])

dnl search for the fastjet-config script
if test "$with_fastjet" = ""; then
   AC_PATH_PROG(fjconfig, fastjet-config, no)
else
   AC_PATH_PROG(fjconfig, fastjet-config, no, ${with_fastjet}/bin)
fi

if test "${fjconfig}" = "no"; then
   AC_MSG_CHECKING(FastJet)
   AC_MSG_RESULT(no);
   $2
else

   dnl now see if FastJet is functional
   save_CXXFLAGS="$CXXFLAGS"
   save_LIBS="$LIBS"

   CXXFLAGS="${CXXFLAGS} `${fjconfig} --cxxflags`"
   LIBS="${LIBS} `${fjconfig} --libs --plugins`"

   AC_MSG_CHECKING([if FastJet is functional])
   AC_LANG_PUSH(C++)
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <fastjet/ClusterSequence.hh>
   ]], [[
fastjet::PseudoJet pj=fastjet::PtYPhiM(10.0,0.5,1.0,0.0);
   ]])], [fjok='yes'], [fjok='no'])
   AC_MSG_RESULT([$fjok])
   AC_LANG_POP()
   CXXFLAGS="$save_CXXFLAGS"
   LIBS="$save_LIBS"

   AC_MSG_CHECKING(FastJet)
   if test "${fjok}" = "yes"; then
      FASTJET_CXXFLAGS="`${fjconfig} --cxxflags`"
      FASTJET_LIBS="`${fjconfig} --libs --plugins`"
      AC_MSG_RESULT(yes)
      $1
   else
      AC_MSG_RESULT(no)
      $2
   fi
fi
])

dnl CHECK FASTJET END
