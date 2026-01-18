dnl gosam.m4 -- checks for gosam package and required helper packages
dnl

AC_DEFUN([WO_PROG_GOSAM],
[dnl
AC_ARG_ENABLE([gosam],
  [AS_HELP_STRING([--enable-gosam],
     [(experimental) enable GoSam for NLO matrix elements [[no]]])],
  [], [enable_gosam="no"])

AC_ARG_WITH([gosam],
  [AS_HELP_STRING([--with-gosam=dir],
     [assume the given directory for GoSam])])

AC_ARG_WITH([golem],
  [AS_HELP_STRING([--with-golem=dir],
     [assume the given directory for Golem])])

AC_ARG_WITH([form],
  [AS_HELP_STRING([--with-form=dir],
     [assume the given directory for Form])])

AC_ARG_WITH([qgraf],
  [AS_HELP_STRING([--with-qgraf=dir],
     [assume the given directory for QGRAF])])

AC_ARG_WITH([ninja],
  [AS_HELP_STRING([--with-ninja=dir],
     [assume the given directory for Ninja])])

AC_ARG_WITH([samurai],
  [AS_HELP_STRING([--with-samurai=dir],
     [assume the given directory for Samurai])])

unset GOSAM_DIR

if test "$enable_gosam" = "yes"; then

  if test "$with_gosam" = ""; then
    AC_PATH_PROG(gosam_exe, [gosam.py], [no])
  else
    AC_PATH_PROG(gosam_exe, [gosam.py], no, ${with_gosam}/bin)
  fi

  if test "$gosam_exe" = "no"; then
    AC_MSG_ERROR([GoSam is enabled but not found])
  else
    gosam_bindir=`dirname $gosam_exe`
    GOSAM_DIR=`dirname $gosam_bindir`
    echo "Gosam dir is " $GOSAM_DIR
  fi

  AC_MSG_CHECKING([the GoSam version])
  wo_gosam_version=`$gosam_exe --version | $GREP "(rev" | $SED 's/GoSam //g' | $SED 's/ (rev.*$//g'`
  GOSAM_VERSION=$wo_gosam_version
  AC_MSG_RESULT([$wo_gosam_version])
  AC_SUBST([GOSAM_VERSION])

  save_path=$PATH
  save_ld_library_path=$LD_LIBRARY_PATH

  AC_MSG_CHECKING([for gosam_setup_env.sh])
  gosam_env=${GOSAM_DIR}/bin/gosam_setup_env.sh
  
  if test -f $gosam_env; then
    AC_MSG_RESULT([$gosam_env])
    . $gosam_env
  else
    AC_MSG_RESULT([no])
    PATH=${GOSAM_DIR}/bin:$PATH
    LD_LIBRARY_PATH=${GOSAM_DIR}/lib:${GOSAM_DIR}/lib64:$LD_LIBRARY_PATH
  fi

  WO_PROG_GOLEM()
  if test -z "$GOLEM_DIR"; then
    AC_MSG_ERROR([GoSam is enabled but Golem is not found])
  fi

  WO_PROG_FORM()
  if test -z "$FORM_DIR"; then
    AC_MSG_ERROR([GoSam is enabled but Form is not found])
  fi

  WO_PROG_QGRAF()
  if test -z "$QGRAF_DIR"; then
    AC_MSG_ERROR([GoSam is enabled but QGRAF is not found])
  fi

  WO_PROG_NINJA()
  if test -z "$NINJA_DIR"; then
    AC_MSG_ERROR([GoSam is enabled but Ninja is not found])
  fi

  WO_PROG_SAMURAI()
  if test -z "$SAMURAI_DIR"; then
    AC_MSG_ERROR([GoSam is enabled but Samurai is not found])
  fi

  PATH=$save_path
  LD_LIBRARY_PATH=$save_ld_library_path

else

  AC_MSG_CHECKING([for GoSam])
  AC_MSG_RESULT([(disabled)])

fi

AC_SUBST([GOSAM_DIR])

if test "$enable_gosam" = "yes"; then
   GOSAM_AVAILABLE_FLAG=".true."
else
   GOSAM_AVAILABLE_FLAG=".false."
fi
AC_SUBST([GOSAM_AVAILABLE_FLAG])

AM_CONDITIONAL([GOSAM_AVAILABLE], [test "$enable_gosam" = "yes"])

]) dnl WO_PROG_GOSAM


AC_DEFUN([WO_PROG_GOLEM],
[dnl
  unset GOLEM_DIR
  if test -n "$with_golem"; then
    echo "Checking for golem in " ${with_golem}/lib
    WO_PATH_LIB(golem_lib, [golem], [libgolem.la], ${with_golem}/lib:${with_golem}/lib64)
  else
    WO_PATH_LIB(golem_lib, [golem], [libgolem.la], $LD_LIBRARY_PATH)
  fi
  if test "$golem_lib" != "no"; then
    golem_libdir=`dirname $golem_lib`
    GOLEM_DIR=`dirname $golem_libdir`
  fi
  AC_SUBST([GOLEM_DIR])
])

AC_DEFUN([WO_PROG_FORM],
[dnl
  unset FORM_DIR
  if test -n "$with_form"; then
    AC_PATH_PROG(form_exe, [form], no, ${with_form}/bin)
  else
    AC_PATH_PROG(form_exe, [form], no)
  fi
  if test "$form_exe" != "no"; then
    form_bindir=`dirname $form_exe`
    FORM_DIR=`dirname $form_bindir`
  fi
  AC_SUBST([FORM_DIR])
])

AC_DEFUN([WO_PROG_QGRAF],
[dnl
  unset QGRAF_DIR
  if test -n "$with_qgraf"; then
    AC_PATH_PROG(qgraf_exe, [qgraf], no, ${with_qgraf}/bin)
  else
    AC_PATH_PROG(qgraf_exe, [qgraf], no)
  fi
  if test "$qgraf_exe" != "no"; then
    qgraf_bindir=`dirname $qgraf_exe`
    QGRAF_DIR=`dirname $qgraf_bindir`
  fi
  AC_SUBST([QGRAF_DIR])
])

AC_DEFUN([WO_PROG_NINJA],
[dnl
  unset NINJA_DIR
  if test -n "$with_ninja"; then
    echo "Checking for ninja in " ${with_ninja}/lib
    WO_PATH_LIB(ninja_lib, [ninja], [libninja.la], ${with_ninja}/lib:${with_ninja}/lib64)
  else
    WO_PATH_LIB(ninja_lib, [ninja], [libninja.la], $LD_LIBRARY_PATH)
  fi
  if test "$ninja_lib" != "no"; then
    ninja_libdir=`dirname $ninja_lib`
    NINJA_DIR=`dirname $ninja_libdir`
  fi
  AC_SUBST([NINJA_DIR])
])

AC_DEFUN([WO_PROG_SAMURAI],
[dnl
  unset SAMURAI_DIR
  if test -n "$with_samurai"; then
    echo "Checking for samurai in " ${with_samurai}/lib
    WO_PATH_LIB(samurai_lib, [samurai], [libsamurai.la], ${with_samurai}/lib:${with_samurai}/lib64)
  else
    WO_PATH_LIB(samurai_lib, [samurai], [libsamurai.la], $LD_LIBRARY_PATH)
  fi
  if test "$samurai_lib" != "no"; then
    samurai_libdir=`dirname $samurai_lib`
    SAMURAI_DIR=`dirname $samurai_libdir`
  fi
  AC_SUBST([SAMURAI_DIR])
])

