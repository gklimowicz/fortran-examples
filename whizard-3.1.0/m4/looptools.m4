dnl looptools.m4 -- checks for LoopTools library
dnl

include([aux.m4])

### Determine paths to LoopTools components
### If successful, set the conditional LOOPTOOLS_AVAILABLE
### Also: LOOPTOOLS_VERSION LOOPTOOLS_INCLUDES LDFLAGS_LOOPTOOLS
AC_DEFUN([WO_PROG_LOOPTOOLS],
[dnl
AC_REQUIRE([AC_PROG_FC])

AC_ARG_ENABLE([looptools],
  [AS_HELP_STRING([--enable-looptools],
    [enable LoopTools loop integral library [[no]]])],
  [], [enable_looptools="no"])

if test "$enable_looptools" = "yes"; then

  # Guessing the most likely paths
  wo_looptools_path="/usr/local/lib64:/usr/lib64:/opt/local/lib64"

  WO_PATH_LIB(LOOPTOOLS, ooptools, libooptools.a, $wo_looptools_path:$LOOPTOOLS_DIR)

  if test "$LOOPTOOLS_DIR" = ""; then
    enable_looptools="no"
  else
    wo_looptools_includes="-I$LOOPTOOLS_DIR/../include"
  fi

  if test "$enable_looptools" = "yes"; then
    wo_looptools_libdir="-L$LOOPTOOLS_DIR"
    AC_LANG_PUSH([Fortran])
    AC_CHECK_LIB([ooptools],[ltini],
      [LDFLAGS_LOOPTOOLS="$wo_looptools_libdir -looptools"],
      [dnl
    AC_MSG_NOTICE([warning:  ********************************************************])
    AC_MSG_NOTICE([warning:  It seems your LoopTools is too old (v2.5 or older) or   ])
    AC_MSG_NOTICE([warning:  was not compiled properly or compiled with a different  ])
    AC_MSG_NOTICE([warning:  FORTRAN compiler and you forgot to add the proper       ])
    AC_MSG_NOTICE([warning:  runtime to LIBS / LD_LIBRARY_PATH:                      ])
    AC_MSG_NOTICE([warning:  Disabling LoopTools support...                          ])
    AC_MSG_NOTICE([warning:  ********************************************************])
    enable_looptools=no],[$wo_looptools_libdir])
    if test "$enable_looptools" = "yes"; then
      LOOPTOOLS_INCLUDES=$wo_looptools_includes      
    fi
    AC_LANG_POP()

  else
    AC_MSG_CHECKING([for LoopTools])
    AC_MSG_RESULT([(disabled)])
  fi

else
  AC_MSG_CHECKING([for LoopTools])
  AC_MSG_RESULT([(disabled)])
fi

AC_SUBST([LOOPTOOLS_INCLUDES])
AC_SUBST([LDFLAGS_LOOPTOOLS])

if test "$enable_looptools" = "yes"; then
  LOOPTOOLS_AVAILABLE_FLAG=".true."
else
  LOOPTOOLS_AVAILABLE_FLAG=".false."
fi
AC_SUBST([LOOPTOOLS_AVAILABLE_FLAG])

AM_CONDITIONAL([LOOPTOOLS_AVAILABLE], [test "$enable_looptools" = "yes"])
])
