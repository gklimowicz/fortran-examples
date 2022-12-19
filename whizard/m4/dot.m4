AC_DEFUN([AC_PROG_DOT],
[dnl
AC_PATH_PROG(DOT,dot,false)
AM_CONDITIONAL([DOT_AVAILABLE],[test "$DOT" != "false"])
if test "$enable_distribution" = "yes"; then
  if test "$DOT" = "false"; then
    AC_MSG_NOTICE([error: ********************************])
    AC_MSG_NOTICE([error: Dot (graphviz) is not installed.])
    AC_MSG_ERROR([********************************])
  else
    AC_MSG_CHECKING([for dot version])  
    DOTVERSION=`$DOT -V 2>&1 | $SED -n -e 's|.*version* *\(.*\) (.*)$|\1|p'`
    AC_MSG_RESULT([$DOTVERSION])
    AC_CACHE_VAL([wo_dot_cv_integer_version],
      [wo_dot_cv_integer_version="`echo "$DOTVERSION" | \
        $AWK 'NR==1 {
          changequote(<<,>>)dnl
            split (<<$>>1, version, "[.+]+");
            printf ("%d%02d%d", version[1], version[2], version[3])}'`"
          changequote([,])])
    DOTINTEGERVERSION=$wo_dot_cv_integer_version
    AC_MSG_CHECKING([for dot version 2400])
    if test $DOTINTEGERVERSION -ge "2400"; then
       AC_MSG_RESULT([ok])
    else
       AC_MSG_RESULT([<= 2.40.0])
       AC_MSG_NOTICE([error: ***********************************])
       AC_MSG_NOTICE([error: found version $DOTVERSION, too old!])
       AC_MSG_ERROR([***********************************])
    fi    
    AC_SUBST([DOTVERSION])
    AC_SUBST([DOTINTEGERVERSION])
  fi
fi
])
