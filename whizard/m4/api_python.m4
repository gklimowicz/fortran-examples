dnl api_python.m4 -- checks for Python/Python API library
dnl

AC_DEFUN([WO_PROG_PYTHON_API],
[dnl
AC_REQUIRE([AX_PYTHON_DEVEL])

AC_ARG_ENABLE([python],
  [AS_HELP_STRING([--enable-python],
    [enable PYTHON/Cython API for WHIZARD [[no]]])],
  [], [enable_python="no"])

AC_PATH_PROG(cython_exe,[cython],[no])
AC_PATH_PROG(cython3_exe,[cython3],[no])

if test "$enable_python" = "yes"; then
  AC_MSG_CHECKING([for PYTHON API])
  if test "$PYTHON_API" = "yes"; then
    if test "$cython3_exe" = "no" -a "$cython_exe" = "no"; then
      enable_python="no"
      AC_MSG_RESULT([(disabled)])
      AC_MSG_WARN([******************
The Python interface requires Cython.
It will be disabled. ****************])
    else
      AC_MSG_RESULT([(enabled)])
    fi
  else
    enable_python="no"
    AC_MSG_RESULT([(disabled)])
    AC_MSG_WARN([**********************
The Python interface requires Python 3.5.
It will be disabled. ********************])
  fi
else
  AC_MSG_CHECKING([for PYTHON API])
  AC_MSG_RESULT([(disabled)])
fi

if test "$enable_python" = "yes"; then
   PYTHON_API_AVAILABLE_FLAG=".true."
else
   PYTHON_API_AVAILABLE_FLAG=".false."
fi
AC_SUBST([PYTHON_API_AVAILABLE_FLAG])

AM_CONDITIONAL([PYTHON_API_AVAILABLE], [test "$enable_python" = "yes"])
])

