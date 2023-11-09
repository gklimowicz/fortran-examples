dnl lcio.m4 -- checks for LCIO library
dnl

### Determine paths to LCIO components
### If successful, set the conditional LCIO_AVAILABLE
### Also: LCIO_VERSION LCIO_INCLUDES LDFLAGS_LCIO
AC_DEFUN([WO_PROG_LCIO],
[dnl
AC_REQUIRE([AC_PROG_CXX])
AC_REQUIRE([AC_PROG_FC])

AC_ARG_ENABLE([lcio],
  [AS_HELP_STRING([--enable-lcio],
    [enable LCIO for handling event data [[yes]]])],
  [], [enable_lcio="yes"])

# First test for LCIO, then for LCIO_DIR; LCIO takes precedence.
if test -n "$LCIO"; then
  wo_lcio_includes="-I$LCIO/include"
elif test -n "$LCIO_DIR"; then
  wo_lcio_includes="-I$LCIO_DIR/include"
fi

if test "$enable_lcio" = "yes"; then
  AC_MSG_CHECKING([the LCIO version])
  AC_LANG([C++])
  wo_cxxflags_tmp=$CXXFLAGS
  CXXFLAGS="$CXXFLAGS $wo_lcio_includes"
  AC_LINK_IFELSE([dnl
    AC_LANG_PROGRAM([[#include "lcio.h"]],
      [[std::cout << LCIO_MAJOR_VERSION << "." << LCIO_MINOR_VERSION << "." << LCIO_PATCH_LEVEL;]])],
    [dnl
    wk_lcio_version=`./conftest`
    AC_MSG_RESULT([$wk_lcio_version])],
    [dnl
    AC_MSG_RESULT([unknown])
    enable_lcio="no"])
  CXXFLAGS=$wo_cxxflags_tmp
fi
LCIO_VERSION=$wk_lcio_version
AC_SUBST([LCIO_VERSION])

wo_lcio_ldflags="-llcio"

if test "$enable_lcio" = "yes"; then
  wo_require_stdcpp="yes"
  AC_MSG_CHECKING([for LCEventImpl class in -llcio])
# First test for LCIO, then for LCIO_DIR; LCIO takes precedence.
  if test -n "$LCIO"; then
    wo_lcio_ldflags="-Wl,-rpath,$LCIO/lib -Wl,-rpath,$LCIO/lib64 -L$LCIO/lib -L$LCIO/lib64 $wo_lcio_ldflags"
  elif test -n "$LCIO_DIR"; then
    wo_lcio_ldflags="-Wl,-rpath,$LCIO_DIR/lib -Wl,-rpath,$LCIO_DIR/lib64 -L$LCIO_DIR/lib -L$LCIO_DIR/lib64 $wo_lcio_ldflags"
  fi
  wo_libs_tmp=$LIBS
  LIBS="$wo_lcio_ldflags $wo_libs_tmp"
  AC_LANG([C++])
  wo_cxxflags_tmp=$CXXFLAGS
  CXXFLAGS="$CXXFLAGS $wo_lcio_includes"
  AC_LINK_IFELSE([dnl
    AC_LANG_PROGRAM([[#include "IMPL/LCEventImpl.h"]],
      [[using namespace IMPL;  LCEventImpl* evt = new LCEventImpl();]])],
    [],
    [enable_lcio="no"])
  AC_MSG_RESULT([$enable_lcio])
  CXXFLAGS=$wo_cxxflags_tmp
  LIBS=$wo_libs_tmp
else
  AC_MSG_CHECKING([for LCIO])
  AC_MSG_RESULT([(disabled)])
fi

if test "$enable_lcio" = "yes"; then
  LCIO_INCLUDES=$wo_lcio_includes
  LDFLAGS_LCIO=$wo_lcio_ldflags
fi

LCIO_AVAILABLE_FLAG=$enable_lcio

AC_SUBST([LCIO_INCLUDES])
AC_SUBST([LDFLAGS_LCIO])
AC_SUBST([LCIO_AVAILABLE_FLAG])

AM_CONDITIONAL([LCIO_AVAILABLE], [test "$enable_lcio" = "yes"])
])
