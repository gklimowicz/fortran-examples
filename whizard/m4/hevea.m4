dnl noweb.m4 -- checks for HEVEA programs
dnl

### Determine paths to hevea components
AC_DEFUN([WO_PROG_HEVEA],
[dnl
AC_PATH_PROG([HEVEA], [hevea])
AC_PATH_PROG([HACHA], [hacha])
AC_PATH_PROG([IMAGEN], [imagen])
AM_CONDITIONAL([HEVEA_AVAILABLE],
  [test -n "$HEVEA" -a -n "$HACHA" -a -n "$IMAGEN"])
])
