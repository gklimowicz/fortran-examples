dnl hoppet.m4 -- checks for HOPPET library
dnl

include([aux.m4])

### Determine paths to HOPPET components
### If successful, set the conditional HOPPET_AVAILABLE
### Also: HOPPET_INCLUDES LDFLAGS_HOPPET
###  relies on HOPPET v1.1.3 or newer 
###      (having a hoppet-config script)
AC_DEFUN([WO_PROG_HOPPET],
[dnl
AC_REQUIRE([AC_PROG_FC])

AC_ARG_ENABLE([hoppet],
  [AS_HELP_STRING([--enable-hoppet],
    [enable HOPPET for b quark pdf matching [[no]]])],
  [], [enable_hoppet="no"])

if test "$enable_hoppet" = "yes"; then
  if test -n "$HOPPET_DIR"; then 
    wo_hoppet_config_path=$HOPPET_DIR/bin:$PATH
  else
    wo_hoppet_config_path=$PATH
  fi 
  AC_PATH_PROG([HOPPET_CONFIG], [hoppet-config], [no], 
    [$wo_hoppet_config_path])
  wo_hoppet_path="/usr/local/lib:/usr/lib:/opt/local/lib"
 
  if test "$HOPPET_CONFIG" != "no"; then
    HOPPET_ROOT=`$HOPPET_CONFIG --prefix`
  else  
    enable_hoppet="no"
  fi

  if test "$enable_hoppet" = "yes"; then
    wo_hoppet_includes="-I$HOPPET_ROOT/include/hoppet"
    wo_hoppet_libdir="-L$HOPPET_ROOT/lib -L$HOPPET_ROOT/lib64"
    AC_LANG_PUSH([Fortran])
    AC_CHECK_LIB([hoppet_v1],[hoppetAssign],
      [LDFLAGS_HOPPET="$wo_hoppet_libdir -lhoppet_v1"],
      [dnl
    AC_MSG_NOTICE([warning:  ********************************************************])
    AC_MSG_NOTICE([warning:  It seems your HOPPET was not compiled properly or       ])
    AC_MSG_NOTICE([warning:  compiled with a different FORTRAN compiler and you      ])
    AC_MSG_NOTICE([warning:  forgot to add the proper runtime to                     ])
    AC_MSG_NOTICE([warning:  LIBS / LD_LIBRARY_PATH. Disabling HOPPET support...     ])
    AC_MSG_NOTICE([warning:  ********************************************************])
    enable_hoppet=no],[$wo_hoppet_libdir])
    if test "$enable_hoppet" = "yes"; then
      HOPPET_INCLUDES=$wo_hoppet_includes      
    fi
    AC_LANG_POP()
  else
    AC_MSG_CHECKING([for HOPPET])
    AC_MSG_RESULT([(disabled)])
  fi

  AC_MSG_CHECKING([the HOPPET version])
  wo_hoppet_version=`$HOPPET_CONFIG 2>&1 | $GREP "This is" | $SED 's/.*with hoppet //g' | $SED 's/\.$//g'`
  HOPPET_VERSION=$wo_hoppet_version
  AC_MSG_RESULT([$wo_hoppet_version])
  AC_SUBST([HOPPET_VERSION])

else
  AC_MSG_CHECKING([for HOPPET])
  AC_MSG_RESULT([(disabled)])
fi

AC_SUBST([HOPPET_INCLUDES])
AC_SUBST([LDFLAGS_HOPPET])

if test "$enable_hoppet" = "yes"; then
  HOPPET_AVAILABLE_FLAG=".true."
else
  HOPPET_AVAILABLE_FLAG=".false."
fi
AC_SUBST([HOPPET_AVAILABLE_FLAG])

AM_CONDITIONAL([HOPPET_AVAILABLE], [test "$enable_hoppet" = "yes"])
])
