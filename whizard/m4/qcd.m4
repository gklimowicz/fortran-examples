dnl qcd.m4 -- checks for qcd setup (shower, PYTHIA, MPI)
dnl

include('aux.m4')

AC_DEFUN([WO_PROG_QCD],
[dnl
AC_REQUIRE([AC_PROG_FC])

### Choice to enable or disable (internal) PYTHIA6 package

AC_ARG_ENABLE([pythia6],
  [AS_HELP_STRING([--enable-pythia6],
    [enable internal PYTHIA6 for hadronization [[yes]]])],
  [], [enable_pythia6="yes"])

AC_CACHE_CHECK([whether we want to enable PYTHIA6], 
[wo_cv_pythia6],
[dnl
if test "$enable_pythia6" = "yes"; then
  wo_cv_pythia6=yes
else
  wo_cv_pythia6=no
fi])

AC_ARG_ENABLE([pythia6_eh],
  [AS_HELP_STRING([--enable-pythia6_eh],
    [PYTHIA6 patches for high-energy eh collisions [[no]]])],
    [], [enable_pythia6_eh="no"])

if test "$enable_pythia6" = "yes"; then
  PYTHIA6_AVAILABLE_FLAG=".true."
  AC_MSG_CHECKING([for PYTHIA6])
  AC_MSG_RESULT([(enabled)])
  AC_MSG_CHECKING([for PYTHIA6 eh settings])
  if test "$enable_pythia6_eh" = "yes"; then
    AC_MSG_RESULT([(enabled)])
    PYTHIA6_EH_AVAILABLE_FLAG="yes"
  else
    AC_MSG_RESULT([(disabled)])
    PYTHIA6_EH_AVAILABLE_FLAG="no"
  fi
else
  PYTHIA6_AVAILABLE_FLAG=".false."
  AC_MSG_CHECKING([for PYTHIA6])
  AC_MSG_RESULT([(disabled)])
  if test "$enable_pythia6_eh" = "yes"; then
    AC_MSG_WARN([please enable PYTHIA6 together with the eh switch])
  fi
  PYTHIA6_EH_AVAILABLE_FLAG="no"
  enable_pythia6_eh="no"
fi
AC_SUBST(PYTHIA6_AVAILABLE_FLAG)
AC_SUBST(PYTHIA6_EH_AVAILABLE_FLAG)

AM_CONDITIONAL([PYTHIA6_AVAILABLE], 
   [test "$PYTHIA6_AVAILABLE_FLAG" = ".true."])
AM_CONDITIONAL([PYTHIA6_IS_EH], [test "$enable_pythia6_eh" = "yes"])
])
