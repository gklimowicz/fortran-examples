dnl stdcpp.m4 -- take care of the C++ standard library
dnl if it is required
dnl

AC_DEFUN([WO_PROG_STDCPP],
[dnl

if test "$wo_require_stdcpp" = "yes"; then

  AC_REQUIRE([AC_PROG_CXX])

  ### Checking for static C++ libraries for the static version
  ### This is only necessary for MAC OS X and BSD-like OS
  case $host in
    *-darwin*)
       case "$XCODE_VERSION" in
         1.*|2.*|3.*)
  	wo_ldflags_stdcpp="-lstdc++-static" ;;
         *)
          wo_ldflags_stdcpp="-lstdc++" ;;
       esac ;;
    *-*-freebsd2*|*-*-freebsd3.0*|*-*-freebsdelf3.0*)
	wo_ldflags_stdcpp="-lstdc++-static" ;;
    *)
       wo_ldflags_stdcpp="-lstdc++" ;;
  esac

  AC_MSG_CHECKING([for LDFLAGS_STATIC: host system is $host_os: static flag])
  AC_MSG_RESULT([$wo_ldflags_stdcpp])

else

  AC_MSG_CHECKING([for LDFLAGS_STATIC:])
  AC_MSG_RESULT([(not needed)])

  unset wo_ldflags_stdcpp

fi

LDFLAGS_STATIC="$wo_ldflags_stdcpp"
AC_SUBST([LDFLAGS_STATIC])

])


#####

AC_DEFUN([_AC_PROG_CXX_V_OUTPUT],  
[AC_REQUIRE([AC_PROG_CXX])dnl
AC_LANG_PUSH(C++)dnl
AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])
ac_save_CXXFLAGS=$CXXFLAGS
CXXFLAGS="$CXXFLAGS m4_default([$1], [$ac_cv_prog_cxx_v])"
(eval echo $as_me:__oline__: \"$ac_link\") >&AS_MESSAGE_LOG_FD
ac_cxx_v_output=`eval $ac_link AS_MESSAGE_LOG_FD>&1 2>&1 | grep -v 'Driving:'`
echo "$ac_cxx_v_output" >&AS_MESSAGE_LOG_FD
CXXFLAGS=$ac_save_CXXFLAGS

rm -rf conftest*
AC_LANG_POP(C++)dnl    
# If we are using xlf then replace all the commas with spaces.    
if echo $ac_cxx_v_output | grep xlfentry >/dev/null 2>&1; then
  ac_cxx_v_output=`echo $ac_cxx_v_output | sed 's/,/ /g'`
fi
# On HP/UX there is a line like: "LPATH is: /foo:/bar:/baz" where
# /foo, /bar, and /baz are search directories for the Fortran linker.
# Here, we change these into -L/foo -L/bar -L/baz (and put it first):
ac_cxx_v_output="`echo $ac_cxx_v_output |
        grep 'LPATH is:' |
        sed 's,.*LPATH is\(: *[[^ ]]*\).*,\1,;s,: */, -L/,g'` $ac_cxx_v_output"

# The Intel C++ compiler's output is rather verbose, we only need the linker information, located at the end.
if echo $ac_cxx_v_output | grep 'mGLOB_options_string' >/dev/null 2>&1; then
  ac_cxx_v_output="`echo $ac_cxx_v_output | sed -n -e 's/.*\( ld \)//p'`"
fi

])# _AC_PROG_CXX_V_OUTPUT

#####

AC_DEFUN([_AC_PROG_CXX_V],
[AC_CACHE_CHECK([how to get verbose linking output from $CXX],
                [ac_cv_prog_cxx_v],
[AC_LANG_ASSERT(C++)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
[ac_cv_prog_cxx_v=
# Try some options frequently used verbose output
for ac_verb in -v -verbose --verbose -V -\#\#\#; do
  _AC_PROG_CXX_V_OUTPUT($ac_verb)
  # look for -l* and *.a constructs in the output
  for ac_arg in $ac_cxx_v_output; do
     case $ac_arg in
        [[\\/]]*.a | ?:[[\\/]]*.a | -[[lLRu]]*)
          ac_cv_prog_cxx_v=$ac_verb
          break 2 ;;
     esac
  done
done
if test -z "$ac_cv_prog_cxx_v"; then
   AC_MSG_WARN([cannot determine how to obtain linking information from $CXX])
fi],
                  [AC_MSG_WARN([compilation failed])])
])])# _AC_PROG_CXX_V

#####

AC_DEFUN([AC_CXX_LIBRARY_LDFLAGS],
[AC_LANG_PUSH(C++)dnl
_AC_PROG_CXX_V  
AC_CACHE_CHECK([for C++ libraries], ac_cv_cxxlibs,
[if test "x$CXXLIBS" != "x"; then
  ac_cv_cxxlibs="$CXXLIBS" # Let the user override the test.
else

_AC_PROG_CXX_V_OUTPUT     

ac_cv_cxxlibs=

# Save positional arguments (if any)  
ac_save_positional="$[@]"

set X $ac_cxx_v_output
while test $[@%:@] != 1; do      
  shift
  ac_arg=$[1]
  case $ac_arg in   
        [[\\/]]*.a | ?:[[\\/]]*.a)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_cxxlibs, ,
              ac_cv_cxxlibs="$ac_cv_cxxlibs $ac_arg")
          ;;
        -bI:*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_cxxlibs, ,
             [_AC_LINKER_OPTION([$ac_arg], ac_cv_cxxlibs)])
          ;;
          # Ignore these flags.
        -lang* | -lcrt0.o | -lc | -lgcc | -libmil | -LANG:=*)
          ;;
        -lcrt1.o | -lcrt1.10.[[1-7]].o | -lcrt2.o |-lcrtbegin.o )
          ;;
        -lkernel32)
          test x"$CYGWIN" != xyes && ac_cv_cxxlibs="$ac_cv_cxxlibs $ac_arg"
          ;;
        -[[LRuY]])      
          # These flags, when seen by themselves, take an argument.
          # We remove the space between option and argument and re-iterate
          # unless we find an empty arg or a new option (starting with -)
          case $[2] in
             "" | -*);;       
             *)
                ac_arg="$ac_arg$[2]"
                shift; shift
                set X $ac_arg "$[@]"
                ;;
          esac
          ;;
        -YP,*)
          for ac_j in `echo $ac_arg | sed -e 's/-YP,/-L/;s/:/ -L/g'`; do
            _AC_LIST_MEMBER_IF($ac_j, $ac_cv_cxxlibs, ,
                               [ac_arg="$ac_arg $ac_j"
                               ac_cv_cxxlibs="$ac_cv_cxxlibs $ac_j"]) 
          done
          ;;
        -[[lLR]]*)
          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_cxxlibs, ,
                             ac_cv_cxxlibs="$ac_cv_cxxlibs $ac_arg")
          ;;
          # Ignore everything else.
  esac
done
# restore positional arguments    
set X $ac_save_positional; shift

# We only consider "LD_RUN_PATH" on Solaris systems.  If this is seen,
# then we insist that the "run path" must be an absolute path (i.e. it
# must begin with a "/").
case `(uname -sr) 2>/dev/null` in
   "SunOS 5"*)
      ac_ld_run_path=`echo $ac_cxx_v_output |
                        sed -n 's,^.*LD_RUN_PATH *= *\(/[[^ ]]*\).*$,-R\1,p'`
      test "x$ac_ld_run_path" != x &&
        _AC_LINKER_OPTION([$ac_ld_run_path], ac_cv_cxxlibs)
      ;;
esac
fi # test "x$CXXLIBS" = "x"
])
CXXLIBS="$ac_cv_cxxlibs"             
AC_SUBST(CXXLIBS)
AC_LANG_POP(C++)dnl
])# AC_CXX_LIBRARY_LDFLAGS
