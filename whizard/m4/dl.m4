dnl dl.m4 -- checks for dynamic-linking library
dnl

### Activate dynamic linking: look for 'dlopen'
AC_DEFUN([WO_PROG_DL],
[dnl
AC_ARG_ENABLE([dl],
  [AS_HELP_STRING([--enable-dl],
    [enable libdl for creating/loading process libraries on-the-fly [[yes]]])],
  [], [enable_dl="yes"])

RTLD_LAZY_VALUE=0
RTLD_NOW_VALUE=0
RTLD_GLOBAL_VALUE=0
RTLD_LOCAL_VALUE=0
AC_SUBST([RTLD_LAZY_VALUE])
AC_SUBST([RTLD_NOW_VALUE])
AC_SUBST([RTLD_GLOBAL_VALUE])
AC_SUBST([RTLD_LOCAL_VALUE])

if test "$enable_dl" = "yes"; then
  AC_LANG_PUSH([C])
  AC_SEARCH_LIBS([dlopen], [dl], [], [enable_dl="no"])
  AC_CHECK_HEADERS([dlfcn.h],[],[enable_dl="no"])
  if test "$enable_dl" = "yes"; then
     AC_MSG_CHECKING([for the values of RTLD_LAZY & friends])
     AC_RUN_IFELSE([AC_LANG_SOURCE([dnl
        #include <stdio.h>
        #include <dlfcn.h>
        
        int main () {
        FILE* handle = fopen ("conftest.out", "w");
           if (!handle) return 1;
           if (!fprintf (handle, "RTLD_LAZY_VALUE=%i\n", RTLD_LAZY)) return 1;
           if (!fprintf (handle, "RTLD_NOW_VALUE=%i\n", RTLD_NOW)) return 1;
           if (!fprintf (handle, "RTLD_GLOBAL_VALUE=%i\n", RTLD_GLOBAL)) return 1;
           if (!fprintf (handle, "RTLD_LOCAL_VALUE=%i\n", RTLD_LOCAL)) return 1;
           if (fclose (handle)) return 1;
        }
      ])],[. ./conftest.out],[enable_dl="no"])
     if test "$enable_dl" = "yes"; then
        AC_MSG_RESULT([success])
     else
        AC_MSG_RESULT([failed])
     fi
  fi
  AC_MSG_CHECKING([for $CC flag to produce position-independent code])
  AC_MSG_RESULT([$lt_prog_compiler_pic_CC])
  CFLAGS_PIC=$lt_prog_compiler_pic_CC
  AC_SUBST([CFLAGS_PIC])
  AC_LANG_POP()
else
  AC_MSG_CHECKING([for dl (dynamic linking)])
  CFLAGS_PIC=''
  AC_SUBST([CFLAGS_PIC])
  AC_MSG_RESULT([(disabled)])
fi

AM_CONDITIONAL([DL_AVAILABLE], [test "$enable_dl" = "yes"])


### For make check of an installed WHIZARD we have to take
### care of library interdependencies on MAC OS X
case $host in
  *-darwin*)
     DYLD_FLAGS="DYLD_LIBRARY_PATH=\$(pwd)/../../src/.libs:\$(pwd)/../../src/main/.libs:\$(pwd)/../../omega/src/.libs:\${DYLD_LIBRARY_PATH}; export DYLD_LIBRARY_PATH" ;; 
  *)
     DYLD_FLAGS="" ;;
esac
AC_SUBST(DYLD_FLAGS)
])

AC_DEFUN([WO_CC_FLAGS],
[dnl
])

