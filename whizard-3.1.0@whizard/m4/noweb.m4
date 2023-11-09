dnl noweb.m4 -- checks for NOWEB programs
dnl

### Determine paths to noweb components
AC_DEFUN([WO_PROG_NOWEB],
[dnl
AC_ARG_ENABLE([noweb],
  [AS_HELP_STRING([--disable-noweb],
    [disable the noweb programs, even if available [[no]]])])
if test "$enable_noweb" != "no"; then
AC_PATH_PROG([NOTANGLE], [notangle])
AC_PATH_PROG([CPIF], [cpif])
AC_PATH_PROG([NOWEAVE], [noweave])		
fi
AC_SUBST([NOTANGLE])
AC_SUBST([NOWEAVE])
AC_SUBST([CPIF])
if test "$enable_distribution" = "yes"; then
if test "$NOTANGLE" = "" -o "$CPIF" = "" -o "$NOWEAVE" = ""; then
AC_MSG_NOTICE([error: **************************************])
AC_MSG_NOTICE([error: Noweb not or not completely installed.])
AC_MSG_ERROR([**************************************])
fi
fi
AC_ARG_ENABLE([noweb-force],
  [AS_HELP_STRING([--disable-noweb-force],
    [force to disable the noweb programs, even if available, which is for distribution testing purposes only. The default never has any effect. [[no]]])])
if test "$enable_noweb_force" = "no"; then
AC_MSG_WARN([**************************************************])
AC_MSG_WARN([ Noweb switched off by force for testing purposes.])
AC_MSG_WARN([**************************************************])
fi
AM_CONDITIONAL([NOWEB_AVAILABLE],
  [test "$enable_noweb" != "no" -a "$enable_noweb_force" != "no" -a -n "$NOTANGLE" -a -n "$CPIF" -a -n "$NOWEAVE"])
])



