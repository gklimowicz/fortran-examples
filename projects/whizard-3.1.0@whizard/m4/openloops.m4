dnl openloops.m4 -- checks for OpenLoops package
dnl


AC_DEFUN([WO_PROG_OPENLOOPS],
[dnl
AC_ARG_ENABLE([openloops],
  [AS_HELP_STRING([--enable-openloops],
     [enable OpenLoops for NLO matrix elements [[no]]])],
  [], [enable_openloops="no"])

AC_ARG_WITH([openloops],
  [AS_HELP_STRING([--with-openloops=dir],
     [assume the given directory for OpenLoops])])

unset OPENLOOPS_DIR

if test "$enable_openloops" = "yes"; then

  unset OPENLOOPS_DIR
  if test -n "$with_openloops"; then
    WO_PATH_LIB(openloops_lib, [openloops], [libopenloops.${SHRLIB_EXT}], ${with_openloops}/lib:${with_openloops}/lib64)
  else
    WO_PATH_LIB(openloops_lib, [openloops], [libopenloops.${SHRLIB_EXT}], $LD_LIBRARY_PATH)
  fi
  if test "$openloops_lib" != "no"; then
    openloops_libdir=`dirname $openloops_lib`
    OPENLOOPS_DIR=`dirname $openloops_libdir`
  fi

else

  AC_MSG_CHECKING([for OpenLoops])
  AC_MSG_RESULT([(disabled)])

fi

AC_SUBST([OPENLOOPS_DIR])

if test -n "$OPENLOOPS_DIR"; then
  wo_openloops_includes="-I$OPENLOOPS_DIR/lib_src/openloops/mod"
  wo_openloops_ldflags="-Wl,-rpath,$OPENLOOPS_DIR/lib -L$OPENLOOPS_DIR/lib -Wl,-rpath,$OPENLOOPS_DIR/lib64 -L$OPENLOOPS_DIR/lib64 -lopenloops"
  wo_openloops_versionfile="$OPENLOOPS_DIR/pyol/config/default.cfg"
fi

if test "$enable_openloops" = "yes" -a "$openloops_lib" != "no"; then
    AC_MSG_CHECKING([for standard OpenLoops processes])
    
    if test -f "$OPENLOOPS_DIR/proclib/libopenloops_ppllj_lt.info" && test -f "$OPENLOOPS_DIR/proclib/libopenloops_eett_lt.info" && test -n "`$GREP 'eexttxg' $OPENLOOPS_DIR/proclib/libopenloops_eett_lt.info`" && test -f "$OPENLOOPS_DIR/proclib/libopenloops_tbw_lt.info"; then
       AC_MSG_RESULT([ OpenLoops processes ppllj/eett/tbw are installed])
       OPENLOOPS_AVAILABLE_FLAG=".true."
    else
       AC_MSG_RESULT([ OpenLoops processes ppllj/eett/tbw are not installed])
       AC_MSG_NOTICE([error: *************************************************************])
       AC_MSG_NOTICE([error: OpenLoops standard process is not installed, please install  ])
       AC_MSG_NOTICE([error:    ppllj, eett and tbw with compile_extra=1                  ])
       AC_MSG_NOTICE([error: *************************************************************])
       OPENLOOPS_AVAILABLE_FLAG=".false."
       enable_openloops="no"
    fi
    OPENLOOPS_INCLUDES=$wo_openloops_includes
    LDFLAGS_OPENLOOPS=$wo_openloops_ldflags
    AC_MSG_CHECKING([the OpenLoops version])
    wo_openloops_version=`$GREP 'release = ' $wo_openloops_versionfile | $SED 's/release = //g'`
    OPENLOOPS_VERSION=$wo_openloops_version
    AC_MSG_RESULT([$wo_openloops_version])
    if test "$wo_openloops_version" = "1.0.0" || test "$wo_openloops_version" = "1.1.0" || test "$wo_openloops_version" = "1.2.0" || test "$wo_openloops_version" = "1.3.0" || test "$wo_openloops_version" = "1.3.1" || test "$wo_openloops_version" = "2.0.0" || test "$wo_openloops_version" = "2.1.0"; then
       AC_MSG_NOTICE([error: *************************************************************])
       AC_MSG_NOTICE([error: OpenLoops version $wo_openloops_version is too old.          ])
       AC_MSG_NOTICE([error:    Please install at least v2.1.1 or the public beta version.])
       AC_MSG_NOTICE([error: *************************************************************])
       OPENLOOPS_AVAILABLE_FLAG=".false."
       enable_openloops="no"
    fi
    AC_SUBST([OPENLOOPS_VERSION])
	
else
   OPENLOOPS_AVAILABLE_FLAG=".false."
   enable_openloops="no"
fi

AC_SUBST([OPENLOOPS_AVAILABLE_FLAG])
AC_SUBST([OPENLOOPS_INCLUDES])
AC_SUBST([LDFLAGS_OPENLOOPS])

AM_CONDITIONAL([OPENLOOPS_AVAILABLE], [test "$enable_openloops" = "yes"])

]) dnl WO_PROG_OPENLOOPS
