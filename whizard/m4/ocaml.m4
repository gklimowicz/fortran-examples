dnl autoconf macros for OCaml
dnl
dnl JR changed obsolete macro into AS_HELP_STRING
dnl JR added check for ocaml binary
dnl JR added check for ocamlcp
dnl JR added check for ocamlweb
dnl JR added routine for lablgtk
dnl JR added check for ocaml version
dnl JR added conditional for CAMLP4
dnl
dnl Copyright © 2009      Richard W.M. Jones
dnl Copyright © 2009      Stefano Zacchiroli
dnl Copyright © 2000-2005 Olivier Andrieu
dnl Copyright © 2000-2005 Jean-Christophe Filliâtre
dnl Copyright © 2000-2005 Georges Mariano
dnl
dnl For documentation, please read the ocaml.m4 man page.

AC_DEFUN([AC_PROG_OCAML],
[dnl
AC_ARG_ENABLE([ocaml],		
  [AS_HELP_STRING([--disable-ocaml],
    [disable the OCaml parts, even if OCaml available [[no]]])])
  if test "$enable_ocaml" != "no"; then
     # checking for ocamlc
     AC_PATH_TOOL([OCAML], [ocaml], [])
     AC_PATH_TOOL([OCAMLC],[ocamlc],[no])
	
     if test "$OCAMLC" != "no"; then
        OCAMLVERSION=`$OCAMLC -v | $SED -n -e 's|.*version* *\(.*\)$|\1|p'`
        #####
        # JR inserted this ocamlintegerversion for version checking
        # [tho] made it rubust for OCaml 4.00
        #####
        AC_CACHE_VAL([wo_ocaml_cv_integer_version],
          [wo_ocaml_cv_integer_version="`echo "$OCAMLVERSION" | \
            $AWK 'NR==1 {
              changequote(<<,>>)dnl
                split (<<$>>1, version, "[.+]+");
                printf ("%d%02d%03d", version[1], version[2], version[3])}'`"
              changequote([,])])
        OCAMLINTEGERVERSION=$wo_ocaml_cv_integer_version
        #####
        AC_MSG_RESULT([OCaml version is $OCAMLVERSION])
        OCAMLLIB=`$OCAMLC -where 2>/dev/null || $OCAMLC -v|tail -1|cut -d ' ' -f 4`
        AC_MSG_RESULT([OCaml library path is $OCAMLLIB])
        	
   
        AC_SUBST([OCAMLVERSION])
        AC_SUBST([OCAMLINTEGERVERSION])
        AC_SUBST([OCAMLLIB])
   
        # checking for ocamlopt
        AC_PATH_TOOL([OCAMLOPT],[ocamlopt],[no])
        OCAMLBEST=byte
        if test "$OCAMLOPT" = "no"; then
   	AC_MSG_WARN([Cannot find ocamlopt; bytecode compilation only.])
        else
   	TMPVERSION=`$OCAMLOPT -v | $SED -n -e 's|.*version* *\(.*\)$|\1|p' `
   	if test "$TMPVERSION" != "$OCAMLVERSION" ; then
   	    AC_MSG_RESULT([versions differs from ocamlc; ocamlopt discarded.])
   	    OCAMLOPT=no
   	else
   	    OCAMLBEST=opt
   	fi
        fi
   
        AC_SUBST([OCAMLBEST])
   
        # checking for ocamlc.opt
        AC_PATH_TOOL([OCAMLCDOTOPT],[ocamlc.opt],[no])
        if test "$OCAMLCDOTOPT" != "no"; then
   	TMPVERSION=`$OCAMLCDOTOPT -v | $SED -n -e 's|.*version* *\(.*\)$|\1|p' `
   	if test "$TMPVERSION" != "$OCAMLVERSION" ; then
   	    AC_MSG_RESULT([versions differs from ocamlc; ocamlc.opt discarded.])
   	else
   	    OCAMLC=$OCAMLCDOTOPT
   	fi
        fi
   
        # checking for ocamlopt.opt
        if test "$OCAMLOPT" != "no" ; then
   	AC_PATH_TOOL([OCAMLOPTDOTOPT],[ocamlopt.opt],[no])
   	if test "$OCAMLOPTDOTOPT" != "no"; then
   	   TMPVERSION=`$OCAMLOPTDOTOPT -v | $SED -n -e 's|.*version* *\(.*\)$|\1|p' `
   	   if test "$TMPVERSION" != "$OCAMLVERSION" ; then
   	      AC_MSG_RESULT([version differs from ocamlc; ocamlopt.opt discarded.])
   	   else
   	      OCAMLOPT=$OCAMLOPTDOTOPT
   	   fi
           fi
        fi
   
        AC_SUBST([OCAMLOPT])
     fi
   
     AC_SUBST([OCAMLC])

     # Allow to use flags set by environment variable OCAMLFLAGS
     AC_SUBST([OCAMLFLAGS],[$OCAMLFLAGS])

     # checking for ocamldep
     AC_PATH_TOOL([OCAMLDEP],[ocamldep],[no])
   
     # checking for ocamlmktop
     AC_PATH_TOOL([OCAMLMKTOP],[ocamlmktop],[no])
   
     # checking for ocamlmklib
     AC_PATH_TOOL([OCAMLMKLIB],[ocamlmklib],[no])
   
     # checking for ocamldoc
     AC_PATH_TOOL([OCAMLDOC],[ocamldoc],[no])
   
     # checking for ocamlbuild
     AC_PATH_TOOL([OCAMLBUILD],[ocamlbuild],[no])
  fi
  AM_CONDITIONAL([OCAML_AVAILABLE],
     [test "$enable_ocaml" != "no"])
])
])


AC_DEFUN([AC_PROG_OCAMLLEX],
[dnl
  # checking for ocamllex
  AC_PATH_TOOL([OCAMLLEX],[ocamllex],[no])
  if test "$OCAMLLEX" != "no"; then
    AC_PATH_TOOL([OCAMLLEXDOTOPT],[ocamllex.opt],[no])
    if test "$OCAMLLEXDOTOPT" != "no"; then
	OCAMLLEX=$OCAMLLEXDOTOPT
    fi
  fi
  AC_SUBST([OCAMLLEX])
])

AC_DEFUN([AC_PROG_OCAMLYACC],
[dnl
  AC_PATH_TOOL([OCAMLYACC],[ocamlyacc],[no])
  AC_SUBST([OCAMLYACC])
])

AC_DEFUN([AC_PROG_OCAMLCP],
[dnl
  AC_PATH_TOOL([OCAMLCP],[ocamlcp],[no])
  AC_SUBST([OCAMLCP])
])

AC_DEFUN([AC_PROG_CAMLP4],
[dnl
  AC_REQUIRE([AC_PROG_OCAML])dnl

  # checking for camlp4
  AC_PATH_TOOL([CAMLP4],[camlp4],[no])
  if test "$CAMLP4" != "no"; then
     TMPVERSION=`$CAMLP4 -v 2>&1| $SED -n -e 's|.*version *\(.*\)$|\1|p'`
     if test "$TMPVERSION" != "$OCAMLVERSION" ; then
	AC_MSG_RESULT([versions differs from ocamlc])
        CAMLP4=no
     fi
  fi
  AC_SUBST([CAMLP4])
  AM_CONDITIONAL([CAMLP4_AVAILABLE],
     [test "$CAMLP4" != "no"])

  # checking for companion tools
  AC_PATH_TOOL([CAMLP4BOOT],[camlp4boot],[no])
  AC_PATH_TOOL([CAMLP4O],[camlp4o],[no])
  AC_PATH_TOOL([CAMLP4OF],[camlp4of],[no])
  AC_PATH_TOOL([CAMLP4OOF],[camlp4oof],[no])
  AC_PATH_TOOL([CAMLP4ORF],[camlp4orf],[no])
  AC_PATH_TOOL([CAMLP4PROF],[camlp4prof],[no])
  AC_PATH_TOOL([CAMLP4R],[camlp4r],[no])
  AC_PATH_TOOL([CAMLP4RF],[camlp4rf],[no])
  AC_SUBST([CAMLP4BOOT])
  AC_SUBST([CAMLP4O])
  AC_SUBST([CAMLP4OF])
  AC_SUBST([CAMLP4OOF])
  AC_SUBST([CAMLP4ORF])
  AC_SUBST([CAMLP4PROF])
  AC_SUBST([CAMLP4R])
  AC_SUBST([CAMLP4RF])
])


AC_DEFUN([AC_PROG_FINDLIB],
[dnl
  AC_REQUIRE([AC_PROG_OCAML])dnl

  # checking for ocamlfind
  AC_PATH_TOOL([OCAMLFIND],[ocamlfind],[no])
  AC_SUBST([OCAMLFIND])
])


dnl Thanks to Jim Meyering for working this next bit out for us.
dnl XXX We should define AS_TR_SH if it's not defined already
dnl (eg. for old autoconf).
AC_DEFUN([AC_CHECK_OCAML_PKG],
[dnl
  AC_REQUIRE([AC_PROG_FINDLIB])dnl

  AC_MSG_CHECKING([for OCaml findlib package $1])

  unset found
  unset pkg
  found=no
  for pkg in $1 $2 ; do
    if $OCAMLFIND query $pkg >/dev/null 2>/dev/null; then
      AC_MSG_RESULT([found])
      AS_TR_SH([OCAML_PKG_$1])=$pkg
      found=yes
      break
    fi
  done
  if test "$found" = "no" ; then
    AC_MSG_RESULT([not found])
    AS_TR_SH([OCAML_PKG_$1])=no
  fi

  AC_SUBST(AS_TR_SH([OCAML_PKG_$1]))
])


AC_DEFUN([AC_CHECK_OCAML_MODULE],
[dnl
  AC_MSG_CHECKING([for OCaml module $2])

  cat > conftest.ml <<EOF
open $3
EOF
  unset found
  for $1 in $$1 $4 ; do
    if $OCAMLC -c -I "$$1" conftest.ml >&5 2>&5 ; then
      found=yes
      break
    fi
  done

  if test "$found" ; then
    AC_MSG_RESULT([$$1])
  else
    AC_MSG_RESULT([not found])
    $1=no
  fi
  AC_SUBST([$1])
])


dnl XXX Cross-compiling
AC_DEFUN([AC_CHECK_OCAML_WORD_SIZE],
[dnl
  AC_MSG_CHECKING([for OCaml compiler word size])
  cat > conftest.ml <<EOF
  print_endline (string_of_int Sys.word_size)
  EOF
  OCAML_WORD_SIZE=`ocaml conftest.ml`
  AC_MSG_RESULT([$OCAML_WORD_SIZE])
  AC_SUBST([OCAML_WORD_SIZE])
])

dnl Check for ocamlweb
AC_DEFUN([AC_PROG_OCAMLWEB],
[dnl
  AC_PATH_TOOL([OCAMLWEB],[ocamlweb],[no])
  AC_SUBST([OCAMLWEB])
  if test "$OCAMLWEB" != "no"; then
     OCAMLWEBVERSION=`$OCAMLWEB --version | $SED -n -e 's|.*version* *\(.*\)$|\1|p'`
     AC_CACHE_VAL([wo_ocamlweb_cv_integer_version],
       [wo_ocamlweb_cv_integer_version="`$OCAMLWEB --version | \
         $AWK 'NR==1 && [$]5 ~ /version/ {
           changequote(<<,>>)dnl
             split (<<$>>6, version, "[.+]+");
             printf ("%d%02d%03d", version[1], version[2], version[3])}'`"
           changequote([,])])
     OCAMLWEBINTEGERVERSION=$wo_ocamlweb_cv_integer_version
     if test "$OCAMLINTEGERVERSION" -ge "$1"; then
          AC_MSG_RESULT([OCamlweb version is $OCAMLWEBVERSION: ok])
     else			  
     	  AC_MSG_RESULT([OCamlweb version is $OCAMLWEBVERSION: too old])
	  OCAMLWEB=no
     fi
  fi	  
  AC_SUBST([OCAMLWEB])
  AC_SUBST([OCAMLWEBVERSION])
  AM_CONDITIONAL([OCAMLWEB_AVAILABLE],[test "$OCAMLWEB" != "no"])
if test "$enable_distribution" = "yes"; then
if test "$OCAMLWEB" = "no"; then
AC_MSG_NOTICE([error: **************************])
AC_MSG_NOTICE([error: Ocamlweb is not installed.])
AC_MSG_ERROR([**************************])
fi
fi
])
 
dnl
dnl --------------------------------------------------------------------
dnl
AC_DEFUN([AC_PROG_OCAML_LABLGTK],
[AC_REQUIRE([AC_PROG_OCAML])
AC_SUBST(LABLGTKDIR)
LABLGTKDIR=$OCAMLLIB
AC_MSG_CHECKING([for OCaml/GTK+ toolkit directory])
if test -f $LABLGTKDIR/lablgtk.cma; then
  AC_MSG_RESULT([$LABLGTKDIR])
else
  LABLGTKDIR=$OCAMLLIB/lablgtk
  if test -f $LABLGTKDIR/lablgtk.cma; then
    AC_MSG_RESULT([$LABLGTKDIR])
  else
    AC_MSG_RESULT([not found])
  fi
fi])

dnl
dnl --------------------------------------------------------------------
dnl
AC_DEFUN([AC_OCAML_VERSION_CHECK],
[AC_REQUIRE([AC_PROG_OCAML])
AC_MSG_CHECKING([for OCaml version $1])
if test $OCAMLINTEGERVERSION -ge "$1"; then
   AC_MSG_RESULT([ok])
else
   AC_MSG_RESULT([<= 4.02.2])
   AC_MSG_NOTICE([error: *************************************])
   AC_MSG_NOTICE([error: found version $OCAMLVERSION, too old!])
   AC_MSG_ERROR([*************************************])
fi])

dnl
dnl --------------------------------------------------------------------
dnl
AC_DEFUN([AC_OCAML_BIGARRAY_MODULE],
[AC_REQUIRE([AC_PROG_OCAML])
AC_MSG_CHECKING([for OCaml bigarray implementation])
if test $OCAMLINTEGERVERSION -ge 408000; then
   OCAML_BIGARRAY_COMPAT=bigarray_module
   OCAML_BIGARRAY_CMA=
   OCAML_BIGARRAY_CMXA=
   AC_MSG_RESULT([module])
elif test $OCAMLINTEGERVERSION -ge 406000; then
   OCAML_BIGARRAY_COMPAT=bigarray_module
   OCAML_BIGARRAY_CMA=bigarray.cma
   OCAML_BIGARRAY_CMXA=bigarray.cmxa
   AC_MSG_RESULT([transitioning])
else
   OCAML_BIGARRAY_COMPAT=bigarray_library
   OCAML_BIGARRAY_CMA=bigarray.cma
   OCAML_BIGARRAY_CMXA=bigarray.cmxa
   AC_MSG_RESULT([library])
fi
AC_SUBST([OCAML_BIGARRAY_COMPAT])
AC_SUBST([OCAML_BIGARRAY_CMA])
AC_SUBST([OCAML_BIGARRAY_CMXA])])

