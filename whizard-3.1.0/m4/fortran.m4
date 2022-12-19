dnl fortran.m4 -- Fortran compiler checks beyond Autoconf built-ins
dnl

dnl The standard Fortran compiler test is AC_PROG_FC.
dnl At the end FC, FCFLAGS and FCFLAGS_f90 are set, if successful.

### Determine vendor and version string.
AC_DEFUN([WO_FC_GET_VENDOR_AND_VERSION],
[dnl
AC_REQUIRE([AC_PROG_FC])

AC_CACHE_CHECK([the compiler ID string],
[wo_cv_fc_id_string],
[dnl
$FC -V >conftest.log 2>&1
$FC -version >>conftest.log 2>&1
$FC --version >>conftest.log 2>&1

wo_fc_grep_GFORTRAN=`$GREP -i 'GNU Fortran' conftest.log | head -1`
wo_fc_grep_G95=`$GREP -i 'g95' conftest.log | $GREP -i 'gcc' | head -1`
wo_fc_grep_NAG=`$GREP 'NAG' conftest.log | head -1`
wo_fc_grep_Intel=`$GREP 'IFORT' conftest.log | head -1`
wo_fc_grep_Sun=`$GREP 'Sun' conftest.log | head -1`
wo_fc_grep_Lahey=`$GREP 'Lahey' conftest.log | head -1`
wo_fc_grep_PGF90=`$GREP 'pgf90' conftest.log | head -1`
wo_fc_grep_PGF95=`$GREP 'pgf95' conftest.log | head -1`
wo_fc_grep_PGFORTRAN=`$GREP 'pgfortran' conftest.log | head -1`
wo_fc_grep_PGHPF=`$GREP 'pghpf' conftest.log | head -1`
wo_fc_grep_FLANG=`$GREP 'clang version' conftest.log | head -1`
wo_fc_grep_NVIDIA=`$GREP 'nvfortran' conftest.log | head -1`
wo_fc_grep_default=`cat conftest.log | head -1`

if test -n "$wo_fc_grep_GFORTRAN"; then
  wo_cv_fc_id_string=$wo_fc_grep_GFORTRAN
elif test -n "$wo_fc_grep_G95"; then
  wo_cv_fc_id_string=$wo_fc_grep_G95
elif test -n "$wo_fc_grep_NAG"; then
  wo_cv_fc_id_string=$wo_fc_grep_NAG
elif test -n "$wo_fc_grep_Intel"; then
  wo_cv_fc_id_string=$wo_fc_grep_Intel
elif test -n "$wo_fc_grep_Sun"; then
  wo_cv_fc_id_string=$wo_fc_grep_Sun
elif test -n "$wo_fc_grep_Lahey"; then
  wo_cv_fc_id_string=$wo_fc_grep_Lahey
elif test -n "$wo_fc_grep_PGF90"; then
  wo_cv_fc_id_string=$wo_fc_grep_PGF90
elif test -n "$wo_fc_grep_PGF95"; then
  wo_cv_fc_id_string=$wo_fc_grep_PGF95
elif test -n "$wo_fc_grep_PGFORTRAN"; then
  wo_cv_fc_id_string=$wo_fc_grep_PGFORTRAN
elif test -n "$wo_fc_grep_PGHPF"; then
  wo_cv_fc_id_string=$wo_fc_grep_PGHPF
elif test -n "$wo_fc_grep_FLANG"; then
  wo_cv_fc_id_string=$wo_fc_grep_FLANG
elif test -n "$wo_fc_grep_NVIDIA"; then
  wo_cv_fc_id_string=$wo_fc_grep_NVIDIA
else
  wo_cv_fc_id_string=$wo_fc_grep_default
fi

rm -f conftest.log
])
FC_ID_STRING="$wo_cv_fc_id_string"
AC_SUBST([FC_ID_STRING])

AC_CACHE_CHECK([the compiler vendor],
[wo_cv_fc_vendor],
[dnl
if test -n "$wo_fc_grep_GFORTRAN"; then
  wo_cv_fc_vendor="gfortran"
elif test -n "$wo_fc_grep_G95"; then
  wo_cv_fc_vendor="g95"
elif test -n "$wo_fc_grep_NAG"; then
  wo_cv_fc_vendor="NAG"
elif test -n "$wo_fc_grep_Intel"; then
  wo_cv_fc_vendor="Intel"
elif test -n "$wo_fc_grep_Sun"; then
  wo_cv_fc_vendor="Sun"
elif test -n "$wo_fc_grep_Lahey"; then
  wo_cv_fc_vendor="Lahey"
elif test -n "$wo_fc_grep_PGF90"; then
  wo_cv_fc_vendor="PGI"
elif test -n "$wo_fc_grep_PGF95"; then
  wo_cv_fc_vendor="PGI"
elif test -n "$wo_fc_grep_PGFORTRAN"; then
  wo_cv_fc_vendor="PGI"
elif test -n "$wo_fc_grep_PGHPF"; then
  wo_cv_fc_vendor="PGI"
elif test -n "$wo_fc_grep_FLANG"; then
  wo_cv_fc_vendor="flang"
elif test -n "$wo_fc_grep_NVIDIA"; then
  wo_cv_fc_vendor="NVIDIA"
else
  wo_cv_fc_vendor="unknown"
fi
])
FC_VENDOR="$wo_cv_fc_vendor"


AC_SUBST([FC_VENDOR])

AM_CONDITIONAL([FC_IS_GFORTRAN],
  [test "$FC_VENDOR" = gfortran])

AM_CONDITIONAL([FC_IS_NAG],
  [test "$FC_VENDOR" = NAG])

AC_CACHE_CHECK([the compiler version],
[wo_cv_fc_version],
[dnl
case $FC_VENDOR in
gfortran)
  wo_cv_fc_version=[`echo $FC_ID_STRING | $SED -e 's/[^0-9]*\([0-9]\{1,2\}\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/'`]
  ;;
g95)
  wo_cv_fc_version=[`echo $FC_ID_STRING | $SED -e 's/.*g95 \([0-9][0-9]*\.[0-9][0-9]*\).*$/\1/'`]
  ;;
NAG)
  wo_cv_fc_version=[`echo $FC_ID_STRING | $SED -e 's/.* Release \([0-9][0-9]*\.[0-9][0-9]*\).*$/\1/'`]
  ;;
Intel)
  wo_cv_fc_version=[`echo $FC_ID_STRING | $SED -e 's/.*\([0-9]\{2\}\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/'`]
  ;;
Sun)
  wo_cv_fc_version=[`echo $FC_ID_STRING | $SED -e 's/.* Fortran 95 \([0-9][0-9]*\.[0-9][0-9]*\) .*/\1/'`]
  ;;
PGI)
  wo_cv_fc_version=[`echo $FC_ID_STRING | $SED -e 's/[a-zA-Z\(\)]//g;s/^[0-9]\{2\}//g;s/32.*\|64.*//g'`]
  ;;
flang)
 wo_cv_fc_version=[`echo $FC_ID_STRING | $SED -e 's/.*clang version \([0-9][0-9]*\.[0-9][0-9]*\).*$/\1/'`]
  ;;
NVIDIA)
  wo_cv_fc_version=[`echo $FC_ID_STRING | $SED -e 's/[a-zA-Z\(\)]//g;s/^[0-9]\{2\}//g;s/32.*\|64.*//g'`]
  ;;
*)
  wo_cv_fc_version="unknown"
  ;;
esac
])
FC_VERSION="$wo_cv_fc_version"
AC_SUBST([FC_VERSION])

### Veto old versions of gfortran 4.5/4.6/4.7/4.8/4.9
if test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.5.0" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.5.1" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.5.2" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.5.3" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.5.4" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.6.0" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.6.1" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.6.2" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.6.3" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.6.4" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.7.0" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.7.1" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.7.2" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.7.3" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.7.4" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.8.0"  || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.8.1"  || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.8.2" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.8.3" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.8.4" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.8.5" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.9.0" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.9.1" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.9.2" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.9.3" || test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "4.9.4"; then
FC_IS_GFORTRAN_4="yes"
  else
FC_IS_GFORTRAN_4="no"
fi
### Veto buggy version of gfortran 6.5
if test "$wo_cv_fc_vendor" = "gfortran" -a "$wo_cv_fc_version" = "6.5.0"; then
FC_IS_GFORTRAN_65="yes"
  else
FC_IS_GFORTRAN_65="no"
fi
AC_SUBST([FC_IS_GFORTRAN_4])
AC_SUBST([FC_IS_GFORTRAN_65])
AC_SUBST([FC_IS_NAG])

### Veto old ifort versions 15.0.0/1/2/3/4/5/6/7 and 16.0.0/1/2/3/4 and 17.0.0/1/2/3/4/5/6/7/8
if test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "15.0.0" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "15.0.1" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "15.0.2" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "15.0.3" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "15.0.4" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "15.0.5" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "15.0.6" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "15.0.7" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "16.0.0" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "16.0.1" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "16.0.2" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "16.0.3" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "16.0.4" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "17.0.0" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "17.0.1" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "17.0.2" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "17.0.2" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "17.0.3" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "17.0.4" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "17.0.5" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "17.0.6" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "17.0.7" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "17.0.7" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "17.0.8" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "18.0.0" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "18.0.1" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "18.0.2" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "18.0.2" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "18.0.3" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "18.0.4" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "18.0.5"; then
FC_IS_IFORT15161718="yes"
  else
FC_IS_IFORT15161718="no"
fi
AC_SUBST([FC_IS_IFORT15161718])


### Catch buggy ifort version 19.0.0/1/2
if test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "19.0.0" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "19.0.1" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "19.0.2"; then
FC_IS_IFORT190012="yes"
  else
FC_IS_IFORT190012="no"
fi
AC_SUBST([FC_IS_IFORT190012])

### Catch buggy ifort version 21.1/1/2
if test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "21.0.0" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "21.1.0" || test "$wo_cv_fc_vendor" = "Intel" -a "$wo_cv_fc_version" = "21.2.0"; then
FC_IS_IFORT21012="yes"
  else
FC_IS_IFORT21012="no"
fi
AC_SUBST([FC_IS_IFORT21012])

AC_CACHE_CHECK([the major version],
[wo_cv_fc_major_version],
[wo_cv_fc_major_version=[`echo $wo_cv_fc_version | $SED -e 's/\([0-9][0-9]*\)\..*/\1/'`]
])
FC_MAJOR_VERSION="$wo_cv_fc_major_version"
AC_SUBST([FC_MAJOR_VERSION])

])
### end WO_FC_GET_VENDOR_AND_VERSION

AC_DEFUN([WO_FC_VETO_GFORTRAN_4],
[dnl
if test "$FC_IS_GFORTRAN_4" = "yes"; then
AC_MSG_NOTICE([error: ****************************************])
AC_MSG_NOTICE([error: gfortran 4.X is too old, please upgrade.])
AC_MSG_ERROR([****************************************])
fi 
])

AC_DEFUN([WO_FC_VETO_GFORTRAN_65],
[dnl
if test "$FC_IS_GFORTRAN_65" = "yes"; then
AC_MSG_NOTICE([error: ******************************************************])
AC_MSG_NOTICE([error: gfortran 6.5 is buggy, please use a different version.])
AC_MSG_ERROR([******************************************************])
fi 
])

AC_DEFUN([WO_FC_VETO_IFORT_15_18],
[dnl
if test "$FC_IS_IFORT15161718" = "yes"; then
AC_MSG_NOTICE([error: ***************************************************************])
AC_MSG_NOTICE([error: ifort version < 19 suffers from severe compiler bugs, disabled.])
AC_MSG_ERROR([***************************************************************])
fi 
])

AC_DEFUN([WO_FC_VETO_IFORT_190012],
[dnl
if test "$FC_IS_IFORT190012" = "yes"; then
AC_MSG_NOTICE([error: *************************************************************])
AC_MSG_NOTICE([error: ifort v19.0.0/1/2 suffer from severe compiler bugs, disabled.])
AC_MSG_ERROR([*************************************************************])
fi 
])

AC_DEFUN([WO_FC_VETO_IFORT_21012],
[dnl
if test "$FC_IS_IFORT21012" = "yes"; then
AC_MSG_NOTICE([error: *****************************************************])
AC_MSG_NOTICE([error: ifort v21.0/1/2 suffer from a compiler bug, disabled.])
AC_MSG_ERROR([*****************************************************])
fi
])

### Determine Fortran flags and file extensions
AC_DEFUN([WO_FC_PARAMETERS],
[dnl
AC_REQUIRE([AC_PROG_FC])
AC_REQUIRE([_LT_COMPILER_PIC])
AC_LANG([Fortran])

AC_MSG_CHECKING([for $FC flags])
AC_MSG_RESULT([$FCFLAGS])

AC_MSG_CHECKING([for $FC flag to produce position-independent code])
AC_MSG_RESULT([$lt_prog_compiler_pic_FC])
FCFLAGS_PIC=$lt_prog_compiler_pic_FC
AC_SUBST([FCFLAGS_PIC])

AC_MSG_CHECKING([for $FC source extension])
AC_MSG_RESULT([$ac_fc_srcext])
FC_SRC_EXT=$ac_fc_srcext
AC_SUBST([FC_SRC_EXT])

AC_MSG_CHECKING([for object file extension])
AC_MSG_RESULT([$ac_objext])
OBJ_EXT=$ac_objext
AC_SUBST([OBJ_EXT])
])
### end WO_FC_PARAMETERS


### Determine runtime libraries
### The standard check is insufficient for some compilers
AC_DEFUN([WO_FC_LIBRARY_LDFLAGS],
[dnl
AC_REQUIRE([AC_PROG_FC])
case $host in
*-darwin*)
  system_darwin="yes"
  OS_IS_DARWIN=".true."
  ;;
*)
  system_darwin="no"
  OS_IS_DARWIN=".false."
  ;;
esac
case $FC_VENDOR in
NAG)
  WO_NAGFOR_LIBRARY_LDFLAGS()
  ;;
Intel)
  case $host in
  *-darwin*)
     fcflags_tmp=$FCFLAGS
     FCFLAGS="-shared-intel $FCFLAGS"
     AC_FC_LIBRARY_LDFLAGS()
     fc_libs_tmp=`echo $FCLIBS | $SED -e 's/\/Applications.*/ /'`
     FCFLAGS=$fcflags_tmp
     FCLIBS=$fc_libs_tmp
     ;;
  *)
     AC_FC_LIBRARY_LDFLAGS()
     ;;
  esac
  ;;
*)
  AC_FC_LIBRARY_LDFLAGS()
  ;;
esac
AC_SUBST([OS_IS_DARWIN])
AM_CONDITIONAL([IS_IFORT_DARWIN],
        [test "$system_darwin" = "yes" -a "$FC_VENDOR" = "Intel"])
])

### Check the NAG Fortran compiler
### Use the '-dryrun' feature and extract the libraries from the link command
### Note that the linker is gcc, not ld
AC_DEFUN([WO_NAGFOR_LIBRARY_LDFLAGS],
[dnl
  AC_CACHE_CHECK([Fortran libraries of $FC],
  [wo_cv_fc_libs],
  [dnl
  AC_REQUIRE([WO_FC_CHECK_OPENMP])
  AC_REQUIRE([WO_FC_SET_OPENMP])
  if test -z "$FCLIBS"; then
    AC_LANG([Fortran])
    AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])
    wo_save_fcflags=$FCFLAGS
    if test "$wo_cv_fc_use_openmp" = "yes"; then
      FCFLAGS="-dryrun $wo_cv_fcflags_openmp" 
    else  
      FCFLAGS="-dryrun"
    fi 
    eval "set x $ac_link"
    echo "set x $ac_link"
    shift
    _AS_ECHO_LOG([$[*]])
    wo_nagfor_output=`eval $ac_link AS_MESSAGE_LOG_FD>&1 2>&1`
    echo "$wo_nagfor_output" >&AS_MESSAGE_LOG_FD
    FCFLAGS=$wo_save_fcflags
    wo_cv_fc_libs=`echo $wo_nagfor_output | $SED -e 's/.* -o conftest \(.*\)$/\1/' | $SED -e 's/conftest.$ac_objext //' | $SED -e 's/\/tmp\/conftest.[[0-9]]*\.o/ /' | $SED -e 's/\/tmp\/[[a-zA-Z0-9]]*\/conftest.[[0-9]]*\.o/ /' | $SED -e 's/conftest.o//'`
  else
    wo_cv_fc_libs=$FCLIBS
  fi
  ])
  FCLIBS=$wo_cv_fc_libs
  AC_SUBST([FCLIBS])
])

### Check for basic F95 features
AC_DEFUN([WO_FC_CHECK_F95],
[dnl
AC_CACHE_CHECK([whether $FC supports Fortran 95 features],
[wo_cv_fc_supports_f95],
[dnl
AC_REQUIRE([AC_PROG_FC])
AC_LANG([Fortran])
AC_COMPILE_IFELSE([dnl
  program conftest
    integer, dimension(2) :: ii
    integer :: i
    type :: foo
       integer, pointer :: bar => null ()
    end type foo
    forall (i = 1:2)  ii(i) = i
  contains
    elemental function f(x)
      real, intent(in) :: x
      real :: f
      f = x
    end function f
    pure function g (x) result (gx)
      real, intent(in) :: x
      real :: gx
      gx = x
    end function g
  end program conftest
  ],
  [wo_cv_fc_supports_f95="yes"],
  [wo_cv_fc_supports_f95="no"])
])
FC_SUPPORTS_F95="$wo_cv_fc_supports_f95"
AC_SUBST([FC_SUPPORTS_F95])
if test "$FC_SUPPORTS_F95" = "no"; then
AC_MSG_NOTICE([error: ******************************************************************])
AC_MSG_NOTICE([error: Fortran compiler is not a genuine F95 compiler, configure aborted.])
AC_MSG_ERROR([******************************************************************])
fi])
### end WO_FC_CHECK_F95
 
### Check for the TR15581 extensions (allocatable subobjects)
AC_DEFUN([WO_FC_CHECK_TR15581],
[AC_CACHE_CHECK([whether $FC supports allocatable subobjects],
[wo_cv_fc_allocatable],
[dnl
AC_REQUIRE([AC_PROG_FC])
AC_LANG([Fortran])
AC_COMPILE_IFELSE([dnl
  program conftest
    type :: foo
       integer, dimension(:), allocatable :: bar
    end type foo
  end program conftest
  ],
  [wo_cv_fc_allocatable="yes"],
  [wo_cv_fc_allocatable="no"])
])
FC_SUPPORTS_ALLOCATABLE="$wo_cv_fc_allocatable"
AC_SUBST([FC_SUPPORTS_ALLOCATABLE])
if test "$FC_SUPPORTS_ALLOCATABLE" = "no"; then
AC_MSG_NOTICE([error: ****************************************************************************])
AC_MSG_NOTICE([error: Fortran compiler does not support allocatable structures, configure aborted.])
AC_MSG_ERROR([****************************************************************************])
fi])
### end WO_FC_CHECK_TR15581


### Check for allocatable scalars
AC_DEFUN([WO_FC_CHECK_ALLOCATABLE_SCALARS],
[AC_CACHE_CHECK([whether $FC supports allocatable scalars],
[wo_cv_fc_allocatable_scalars],
[dnl
AC_REQUIRE([AC_PROG_FC])
AC_LANG([Fortran])
AC_COMPILE_IFELSE([dnl
  program conftest
    type :: foo
       integer, allocatable :: bar
    end type foo
  end program conftest
  ],
  [wo_cv_fc_allocatable_scalars="yes"],
  [wo_cv_fc_allocatable_scalars="no"])
])
FC_SUPPORTS_ALLOCATABLE_SCALARS="$wo_cv_fc_allocatable"
AC_SUBST([FC_SUPPORTS_ALLOCATABLE_SCALARS])
])
### end WO_FC_CHECK_ALLOCATABLE_SCALARS


### Check for the C bindings extensions of Fortran 2003
AC_DEFUN([WO_FC_CHECK_C_BINDING],
[AC_CACHE_CHECK([whether $FC supports ISO C binding and standard numeric types],
[wo_cv_fc_c_binding],
[dnl
AC_REQUIRE([AC_PROG_FC])
AC_LANG([Fortran])
AC_COMPILE_IFELSE([dnl
  program conftest
    use iso_c_binding
    type, bind(c) :: t
       integer(c_int) :: i
       real(c_float) :: x1
       real(c_double) :: x2
       complex(c_float_complex) :: z1
       complex(c_double_complex) :: z2
    end type t
  end program conftest
  ],
  [wo_cv_fc_c_binding="yes"],
  [wo_cv_fc_c_binding="no"])
])
FC_SUPPORTS_C_BINDING="$wo_cv_fc_c_binding"
AC_SUBST([FC_SUPPORTS_C_BINDING])
if test "$FC_SUPPORTS_C_BINDING" = "no"; then
AC_MSG_NOTICE([error: *******************************************************************])
AC_MSG_NOTICE([error: Fortran compiler does not support ISO C binding, configure aborted.])
AC_MSG_ERROR([********************************************************************])
fi
])
### end WO_FC_CHECK_C_BINDING


### Check for procedure pointers
AC_DEFUN([WO_FC_CHECK_PROCEDURE_POINTERS],
[AC_CACHE_CHECK([whether $FC supports procedure pointers (F2003)],
[wo_cv_prog_f03_procedure_pointers],
[dnl
AC_REQUIRE([AC_PROG_FC])
AC_LANG([Fortran])
AC_COMPILE_IFELSE([dnl
  program conftest
    type :: foo
      procedure (proc_template), nopass, pointer :: proc => null ()
    end type foo
    abstract interface
      subroutine proc_template ()
      end subroutine proc_template
    end interface
  end program conftest
  ],
  [wo_cv_prog_f03_procedure_pointers="yes"],
  [wo_cv_prog_f03_procedure_pointers="no"])
])
FC_SUPPORTS_PROCEDURE_POINTERS="$wo_cv_prog_f03_procedure_pointers"
AC_SUBST([FC_SUPPORTS_PROCEDURE_POINTERS])
if test "$FC_SUPPORTS_PROCEDURE_POINTERS" = "no"; then
AC_MSG_NOTICE([error: ***************************************************************************])
AC_MSG_NOTICE([error: Fortran compiler does not understand procedure pointers, configure aborted.])
AC_MSG_ERROR([***************************************************************************])
fi])
### end WO_FC_CHECK_PROCEDURE_POINTERS


### Check for the OO extensions of Fortran 2003
AC_DEFUN([WO_FC_CHECK_OO_FEATURES],
[AC_CACHE_CHECK([whether $FC supports OO features (F2003)],
[wo_cv_prog_f03_oo_features],
[dnl
AC_REQUIRE([AC_PROG_FC])
AC_LANG([Fortran])
AC_COMPILE_IFELSE([dnl
  module conftest
    type, abstract :: foo
    contains
       procedure (proc_template), deferred :: proc
    end type foo
    type, extends (foo) :: foobar
    contains
       procedure :: proc
    end type foobar
    abstract interface
      subroutine proc_template (f)
        import foo
        class(foo), intent(inout) :: f
      end subroutine proc_template
    end interface
  contains
    subroutine proc (f)
      class(foobar), intent(inout) :: f
    end subroutine proc
  end module conftest
  program main
    use conftest
  end program main
  ],
  [wo_cv_prog_f03_oo_features="yes"],
  [wo_cv_prog_f03_oo_features="no"])
])
FC_SUPPORTS_OO_FEATURES="$wo_cv_prog_f03_oo_features"
AC_SUBST([FC_SUPPORTS_OO_FEATURES])
])
### end WO_FC_CHECK_OO_FEATURES


### Check for the command line interface of Fortran 2003
### We actually have to link in order to check availability
AC_DEFUN([WO_FC_CHECK_CMDLINE],
[AC_CACHE_CHECK([whether $FC interfaces the command line (F2003)],
  [wo_cv_fc_cmdline],
  [dnl
AC_REQUIRE([AC_PROG_FC])
AC_LANG([Fortran])
AC_LINK_IFELSE([dnl
  program conftest
    call get_command_argument (command_argument_count ())
  end program conftest
  ],
  [wo_cv_fc_cmdline="yes"],
  [wo_cv_fc_cmdline="no"])
])
FC_SUPPORTS_CMDLINE="$wo_cv_fc_cmdline"
AC_SUBST([FC_SUPPORTS_CMDLINE])
if test "$FC_SUPPORTS_CMDLINE" = "no"; then
AC_MSG_NOTICE([error: ******************************************************************])
AC_MSG_NOTICE([error: Fortran compiler does not support get_command_argument; configure aborted.])
AC_MSG_ERROR([******************************************************************])
fi
])
### end WO_FC_CHECK_CMDLINE

### Check whether we can access environment variables (2003 standard)
### We actually have to link in order to check availability
AC_DEFUN([WO_FC_CHECK_ENVVAR],
[AC_CACHE_CHECK([whether $FC provides access to environment variables (F2003)],
  [wo_cv_fc_envvar],
  [dnl
AC_REQUIRE([AC_PROG_FC])
AC_LANG([Fortran])
AC_LINK_IFELSE([dnl
  program conftest
    character(len=256) :: home
    call get_environment_variable ("HOME", home)
  end program conftest
  ],
  [wo_cv_fc_envvar="yes"],
  [wo_cv_fc_envvar="no"])
])
FC_SUPPORTS_ENVVAR="$wo_cv_fc_envvar"
AC_SUBST([FC_SUPPORTS_ENVVAR])
if test "$FC_SUPPORTS_ENVVAR" = "no"; then
AC_MSG_NOTICE([error: ***************************************************************************])
AC_MSG_NOTICE([error: Fortran compiler does not support get_environment_variable; configure aborted.])
AC_MSG_ERROR([***************************************************************************])
fi])
### end WO_FC_CHECK_ENVVAR

### Check whether the flush statement is supported
AC_DEFUN([WO_FC_CHECK_FLUSH],
[AC_CACHE_CHECK([whether $FC supports the flush statement (F2003)],
  [wo_cv_fc_flush],
  [dnl
AC_REQUIRE([AC_PROG_FC])
AC_LANG([Fortran])
AC_COMPILE_IFELSE([dnl
  program conftest
  implicit none
    integer, parameter :: u=50
    open (u, action="readwrite", status="scratch")
    write (u, *) "test"
    flush (u)
    close (u)
  end program conftest
  ],
  [wo_cv_fc_flush="yes"],
  [wo_cv_fc_flush="no"])
])
FC_SUPPORTS_FLUSH="$wo_cv_fc_flush"
AC_SUBST([FC_SUPPORTS_FLUSH])
if test "$FC_SUPPORTS_FLUSH" = "no"; then
AC_MSG_NOTICE([error: ***************************************************************************])
AC_MSG_NOTICE([error: Fortran compiler does not support the flush statement; configure aborted.])
AC_MSG_ERROR([***************************************************************************])
fi])
### end WO_FC_CHECK_FLUSH

### Check for iso_fortran_env
AC_DEFUN([WO_FC_CHECK_ISO_FORTRAN_ENV],
  [AC_CACHE_CHECK([whether $FC supports iso_fortran_env (F2003)],
    [wo_cv_fc_iso_fortran_env],
     AC_REQUIRE([AC_PROG_FC])
     AC_LANG([Fortran])
     [AC_LINK_IFELSE(
        [dnl
        program conftest
        use iso_fortran_env
        implicit none
          integer :: i
          i = input_unit
          i = output_unit
          i = error_unit
          i = iostat_end
          i = iostat_eor
        end program conftest
        ], [wo_cv_fc_iso_fortran_env=yes], [wo_cv_fc_iso_fortran_env=no]
      )]
    )
   if test "$wo_cv_fc_iso_fortran_env" = "no"; then
     AC_MSG_NOTICE([error: *****************************************************************************])
     AC_MSG_NOTICE([error: Fortran compiler does not support iso_fortran_env (F2003); configure aborted.])
     AC_MSG_ERROR([*****************************************************************************])
   fi
  ]
)
### end WO_FC_CHECK_ISO_FORTRAN_ENV

### Check for the TR19767 extensions (submodules)
AC_DEFUN([WO_FC_CHECK_TR19767],
[AC_CACHE_CHECK([whether $FC supports submodules (F2008)],
[wo_cv_fc_submodules],
[dnl
AC_REQUIRE([AC_PROG_FC])
AC_LANG([Fortran])
AC_COMPILE_IFELSE([dnl
  module conftest_mod
    type :: point
       real :: x, y
    end type point
  
    interface
       module function point_dist(a, b) result(distance)
         type(point), intent(in) :: a, b
         real :: distance
       end function point_dist
    end interface
  end module conftest_mod
  
  submodule (conftest_mod) conftest_mod_a
  contains
    module function point_dist(a, b) result(distance)
      type(point), intent(in) :: a, b
      real :: distance
      distance = sqrt((a%x - b%x)**2 + (a%y - b%y)**2)
    end function point_dist
  end submodule conftest_mod_a

  program conftest
    use conftest_mod
  end program conftest
  ],
  [wo_cv_fc_submodules="yes"],
  [wo_cv_fc_submodules="no"])
])
rm -f conftest_mod.*
FC_SUPPORTS_SUBMODULES="$wo_cv_fc_submodules"
AC_SUBST([FC_SUPPORTS_SUBMODULES])
AM_CONDITIONAL([FC_SUBMODULES],
	[test "$wo_cv_fc_submodules" = "yes"])
if test "$FC_SUPPORTS_SUBMODULES" = "no"; then
AC_MSG_NOTICE([NOTE: Fortran compiler does not support submodules (F2008).])
fi])
### end WO_FC_CHECK_TR19767

### Check for wrapping of linker flags 
### (nagfor 'feature': must be wrapped twice)
AC_DEFUN([WO_FC_CHECK_LDFLAGS_WRAPPING],
[AC_CACHE_CHECK([for wrapping of linker flags via -Wl],
  [wo_cv_fc_ldflags_wrapping],
  [dnl
AC_REQUIRE([AC_PROG_FC])
AC_LANG([Fortran])
ldflags_tmp=$LDFLAGS
LDFLAGS=-Wl,-rpath,/usr/lib
AC_LINK_IFELSE(AC_LANG_PROGRAM(),
  [wo_cv_fc_ldflags_wrapping="once"],
  [wo_cv_fc_ldflags_wrapping="unknown"])
if test "$wo_cv_fc_ldflags_wrapping" = "unknown"; then
  LDFLAGS=-Wl,-Wl,,-rpath,,/usr/lib  
  AC_LINK_IFELSE(AC_LANG_PROGRAM(),
    [wo_cv_fc_ldflags_wrapping="twice"])
fi
LDFLAGS=$ldflags_tmp
])
FC_LDFLAGS_WRAPPING="$wo_cv_fc_ldflags_wrapping"
AC_SUBST([FC_LDFLAGS_WRAPPING])
])
### end WO_FC_CHECK_LDFLAGS_WRAPPING

### Check for OpenMP support
AC_DEFUN([WO_FC_CHECK_OPENMP],
[AC_CACHE_CHECK([whether $FC supports OpenMP],
  [wo_cv_fc_openmp],
  [dnl
AC_REQUIRE([AC_PROG_FC])
AC_LANG([Fortran])
case $FC_VENDOR in
gfortran)
  wo_cv_fc_openmp="yes"
  wo_cv_fcflags_openmp="-fopenmp"
  wo_cv_fc_openmp_header="use, intrinsic :: omp_lib"
  ;;
NAG)
  wo_cv_fc_openmp="yes"
  wo_cv_fcflags_openmp="-openmp"
  wo_cv_fc_openmp_header="use, intrinsic :: omp_lib"
  ;;
Intel)
  wo_cv_fc_openmp="yes"
  wo_cv_fcflags_openmp="-qopenmp"
  wo_cv_fc_openmp_header="use :: omp_lib !NODEP!"
  ;;
PGI)
  wo_cv_fc_openmp="yes"
  wo_cv_fcflags_openmp="-mp"
  wo_cv_fc_openmp_header=""
  ;;
NVIDIA)
  wo_cv_fc_openmp="yes"
  wo_cv_fcflags_openmp="-mp"
  wo_cv_fc_openmp_header=""
  ;;
*)
  wo_cv_fc_openmp="no"
  ;;
esac
if test "$wo_cv_fc_openmp"="yes"; then
fcflags_tmp=$FCFLAGS
FCFLAGS="$wo_cv_fcflags_openmp $FCFLAGS"
unset OMP_NUM_THREADS
AC_RUN_IFELSE([dnl
  program conftest
  $wo_cv_fc_openmp_header
  open (11, file="conftest.out", action="write", status="replace")
  write (11, "(I0)") omp_get_max_threads ()
  close (11)
  end program conftest
  ],
  [dnl
  wo_cv_fc_openmp="yes"
  wo_cv_fc_openmp_thread_limit=`cat conftest.out`],
  [wo_cv_fc_openmp="no"],
  [wo_cv_fc_openmp="maybe [cross-compiling]"])
FCFLAGS=$fcflags_tmp
fi
])
if test "$wo_cv_fc_openmp" = "yes"; then
AC_CACHE_CHECK([the default number of threads used by OpenMP],
[wo_cv_fc_openmp_thread_limit])
fi
FC_SUPPORTS_OPENMP="$wo_cv_fc_openmp"
AC_SUBST([FC_SUPPORTS_OPENMP])
])
### end WO_FC_CHECK_OPENMP

### Turn on/off master switch for debugging features
AC_DEFUN([WO_FC_SET_DEBUG],
[dnl
AC_ARG_ENABLE([fc_debug],
  [AS_HELP_STRING([--enable-fc-debug],
    [enable debugging features for the Fortran code [[no]]])])
AC_CACHE_CHECK([whether debugging facilities are enabled], [wo_cv_fc_debug_on],
[dnl
if test "$enable_fc_debug" = "yes"; then
  wo_cv_fc_debug_on=yes
else
  wo_cv_fc_debug_on=no
fi
])
if test "$wo_cv_fc_debug_on" = "yes"; then
  FC_DEBUG_ON=.true.
else
  FC_DEBUG_ON=.false.
fi
AC_SUBST([FC_DEBUG_ON])
])

### Enable/disable OpenMP support
AC_DEFUN([WO_FC_SET_OPENMP],
[dnl
AC_REQUIRE([WO_FC_CHECK_OPENMP])
AC_ARG_ENABLE([fc_openmp],
  [AS_HELP_STRING([--enable-fc-openmp],
    [use OpenMP for the Fortran code [[no]]])])
AC_CACHE_CHECK([whether OpenMP is activated], [wo_cv_fc_use_openmp],
[dnl
if test "$FC_SUPPORTS_OPENMP" = "yes" -a "$enable_fc_openmp" = "yes"; then
  wo_cv_fc_use_openmp="yes"
else
  wo_cv_fc_use_openmp="no"
fi])
AM_CONDITIONAL([FC_USE_OPENMP],
        [test "$wo_cv_fc_use_openmp" = "yes"])
AM_COND_IF([FC_USE_OPENMP],
[FC_OPENMP_ON=""
FC_OPENMP_OFF="!"
FCFLAGS_OPENMP="$wo_cv_fcflags_openmp"
FC_OPENMP_HEADER="$wo_cv_fc_openmp_header"
FC_OPENMP_DEFAULT_MAX_THREADS="$wo_cv_fc_openmp_thread_limit"
],
[FC_OPENMP_ON="!"
FC_OPENMP_OFF=""
FC_OPENMP_DEFAULT_MAX_THREADS="1"
])
AC_SUBST([FC_OPENMP_ON])
AC_SUBST([FC_OPENMP_OFF])
AC_SUBST([FCFLAGS_OPENMP])
AC_SUBST([FC_OPENMP_HEADER])
AC_SUBST([FC_OPENMP_DEFAULT_MAX_THREADS])
])
### end WO_FC_SET_OPENMP

### Enable/disable MPI support (either OpenMPI, MPICH or Intel MPI)
AC_DEFUN([WO_FC_SET_MPI],
[dnl
AC_REQUIRE([WO_FC_FILENAME_CASE_CONVERSION])
AC_ARG_ENABLE([fc_mpi],
  [AS_HELP_STRING([--enable-fc-mpi],
    [use OpenMPI/MPICH/Intel for the Fortran code (default is OpenMPI) [[no]]])],
  [], [enable_fc_mpi="no"])
if test "x$enable_fc_mpi" = "xyes"; then
   if !(test "$FC" == "mpifort" || test "$F77" == "mpifort") \
	&& !(test "$FC" == "mpiifort" || test "$F77" == "mpiifort"); then
        WO_FC_MSG_ERROR_BOX([For MPI please use mpifort or mpiifort (for Intel) as Fortran and F77 compiler!])
   fi
   if test -n "$MPI_DIR"; then
     wo_mpi_config_path=$MPI_DIR/bin:$PATH
   else
     wo_mpi_config_path=$PATH
   fi
   AC_MSG_CHECKING([the requested MPI library])
   AC_ARG_WITH([mpi-lib],
     [  --with-mpi-lib=mpich|openmpi|intel   request an external MPI library.],
     [case "x$withval" in
        x | xno | xyes ) wo_cv_fc_requested_mpilib=openmpi ;;
        * )              wo_cv_fc_requested_mpilib="`echo $withval | $LOWERCASE`" ;;
      esac],
      [wo_cv_fc_requested_mpilib=openmpi])
   case "$wo_cv_fc_requested_mpilib" in
      mpich | openmpi | intel)
        AC_MSG_RESULT([$wo_cv_fc_requested_mpilib])
        ;;
      *)
        AC_MSG_RESULT()
        WO_FC_MSG_ERROR_BOX([argument of --with-mpi-library is $wo_cv_fc_mpilib, but must be one of mpich, openmpi or intel!])
        ;;
   esac
   case "$wo_cv_fc_requested_mpilib" in
      intel )
	AC_PATH_PROG([MPIIFORT],[mpiifort],[no],[$wo_mpi_config_path])
	if test "$MPIIFORT" != "no"; then
	  AC_CACHE_CHECK([the Intel MPI version],
	     [wo_cv_mpiifort_version],
	     [dnl
		wo_cv_mpiifort_version=[`$MPIIFORT --version | head -n 1 | $SED -e 's/ifort (IFORT) \(.*\)\s.*/\1/'`]
		wo_cv_mpiifort_major_version=[`echo $wo_cv_mpiifort_version | $SED -e 's/\([0-9]\)\..*/\1/'`]
	  ])
	  MPI_VERSION=$wo_cv_mpiifort_version
          FCFLAGS_MPI="-lmpifort -lmpi"
	else
	  enable_fc_mpi="no"
        fi
	;;
      mpich )
        AC_PATH_PROG([MPICHVERSION],[mpichversion],[no],[$wo_mpi_config_path])
	if test "$MPICHVERSION" != "no"; then
	  AC_CACHE_CHECK([the MPICH version],
	     [wo_cv_mpich_version],
	     [dnl
	        wo_cv_mpich_version=[`$MPICHVERSION --version | $SED -e 's/MPICH VERSION://'`]
		wo_cv_mpich_major_version=[`echo $wo_cv_mpich_version | $SED -e 's/\([0-9][0-9]*\)\..*/\1/'`]
	  ])
	  MPI_VERSION=$wo_cv_mpich_version
	  FCFLAGS_MPI="-lmpifort"
	else
	  enable_fc_mpi="no"
	fi
	;;
      openmpi )
        AC_PATH_PROG([OMPI_INFO],[ompi_info],[no],[$wo_mpi_config_path])
	if test "$OMPI_INFO" != "no"; then
	  AC_CACHE_CHECK([the OPENMPI version],
	     [wo_cv_openmpi_version],
	     [dnl
	        wo_cv_openmpi_version=[`$OMPI_INFO --version | head -1 | $SED -e 's/Open MPI v//'`]
		wo_cv_openmpi_major_version=[`echo $wo_cv_openmpi_version | $SED -e 's/\([0-9][0-9]*\)\..*/\1/'`]
	  ])
	  MPI_VERSION=$wo_cv_openmpi_version
	  FCFLAGS_MPI="-lmpi -lmpi_usempif08"
	else
	  enable_fc_mpi="no"
	fi
	;;
   esac
else
  AC_MSG_NOTICE([no MPI support demanded])  
fi
MPI_AVAILABLE=$enable_fc_mpi
MPI_LIBRARY=$wo_cv_fc_requested_mpilib
AM_CONDITIONAL([FC_USE_MPI],
        [test "x$enable_fc_mpi" = "xyes"])
])
AC_SUBST([MPI_AVAILABLE])
AC_SUBST([MPI_LIBRARY])
AC_SUBST([MPI_VERSION])
AC_SUBST([FCFLAGS_MPI])
### end WO_FC_SET_MPI


### Check for profiling support
AC_DEFUN([WO_FC_CHECK_PROFILING],
[AC_CACHE_CHECK([whether $FC supports profiling via -pg],
  [wo_cv_fc_profiling],
  [dnl
AC_REQUIRE([AC_PROG_FC])
AC_LANG([Fortran])
fcflags_tmp=$FCFLAGS
FCFLAGS="-pg $FCFLAGS"
rm -f gmon.out
AC_RUN_IFELSE([dnl
  program conftest
  end program conftest
  ],
  [dnl
  if test -f gmon.out; then
    wo_cv_fc_profiling="yes"
  else
    wo_cv_fc_profiling="no"
  fi],
  [wo_cv_fc_profiling="no"],
  [wo_cv_fc_profiling="maybe [cross-compiling]"])
rm -f gmon.out
FCFLAGS=$fcflags_tmp
])
FC_SUPPORTS_PROFILING="$wo_cv_fc_profiling"
AC_SUBST([FC_SUPPORTS_PROFILING])
])
### end WO_FC_CHECK_PROFILING

### Enable/disable profiling support
AC_DEFUN([WO_FC_SET_PROFILING],
[dnl
AC_REQUIRE([WO_FC_CHECK_PROFILING])
AC_ARG_ENABLE([fc_profiling],
  [AS_HELP_STRING([--enable-fc-profiling],
    [use profiling for the Fortran code [[no]]])])
AC_CACHE_CHECK([whether profiling is activated], [wo_cv_fc_prof],
[dnl
if test "$FC_SUPPORTS_PROFILING" = "yes" -a "$enable_fc_profiling" = "yes"; then
  wo_cv_fc_prof="yes"
  FCFLAGS_PROFILING="-pg"       
else
  wo_cv_fc_prof="no"
  FCFLAGS_PROFILING=""
fi])
AC_SUBST(FCFLAGS_PROFILING)
AM_CONDITIONAL([FC_USE_PROFILING],
        [test -n "$FCFLAGS_PROFILING"])
])
### end WO_FC_SET_PROFILING

### Enable/disable impure Omega compilation
AC_DEFUN([WO_FC_SET_OMEGA_IMPURE],
[dnl
AC_REQUIRE([WO_FC_CHECK_F95])
AC_ARG_ENABLE([impure_omega],
  [AS_HELP_STRING([--enable-fc-impure],
    [compile Omega libraries impure [[no]]])])
AC_CACHE_CHECK([the default setting for impure omegalib], [wo_cv_fc_impure],
[dnl
if test "$impure_omega" = "yes" -o "$FC_SUPPORTS_F95" = "no"; then
  wo_cv_fc_impure="yes"
else 
  wo_cv_fc_impure="no"
fi])
AM_CONDITIONAL([FC_IMPURE],
     [test "$wo_cv_fc_impure" = "yes"])
])
### end WO_FC_OMEGA_IMPURE

########################################################################
### Configure kinds.f90
########################################################################

dnl#  splashy error reporting
AC_DEFUN([WO_FC_MSG_ERROR_BOX],
[dnl
  AC_MSG_NOTICE([error: ***************************************************************************])
  AC_MSG_NOTICE([error: $1])
  AC_MSG_ERROR([***************************************************************************])])

AC_DEFUN([WO_FC_MSG_ERROR_BOX2],
[dnl
  AC_MSG_NOTICE([error: ***************************************************************************])
  AC_MSG_NOTICE([error: $1])
  AC_MSG_NOTICE([error: $2])
  AC_MSG_ERROR([***************************************************************************])])

AC_DEFUN([WO_FC_MSG_WARN_BOX],
[dnl
  AC_MSG_WARN([***************************************************************************])
  AC_MSG_WARN([$1])
  AC_MSG_WARN([***************************************************************************])])

dnl#  Check for iso_fortran_env in the incarnation of 2008
AC_DEFUN([WO_FC_CHECK_ISO_FORTRAN_ENV_2008],
 [AC_CACHE_CHECK([whether $FC supports iso_fortran_env (F2008)],
   [wo_cv_fc_iso_fortran_env_2008],
   AC_REQUIRE([AC_PROG_FC])
   AC_LANG([Fortran])
   [AC_LINK_IFELSE(
      [dnl
       program conftest
         use iso_fortran_env
         implicit none
         integer :: i
         i = real32
         i = real64
         i = real128
         i = int8
         i = int16
         i = int32
         i = int64
       end program conftest],
      [wo_cv_fc_iso_fortran_env_2008=yes],
      [wo_cv_fc_iso_fortran_env_2008=no])]
      )
    if test "$wo_cv_fc_iso_fortran_env_2008" = "no"; then
      AC_MSG_NOTICE([error: *****************************************************************************])
      AC_MSG_NOTICE([error: Fortran compiler does not support iso_fortran_env (F2008); configure aborted.])
      AC_MSG_ERROR([*****************************************************************************])
    fi
   ]
)
### end WO_FC_CHECK_ISO_FORTRAN_ENV_2008

dnl#  An extension of ISO_C_BINDING adding definitions of
dnl#  gfortran extensions, flagged as unavailable.
AC_DEFUN([WO_FC_ISO_C_BINDING_GFORTRAN_DUMMY],
 [module iso_c_binding_gfortran
    public
    integer, parameter :: c_float128 = -1
  end module iso_c_binding_gfortran])

dnl#  An empty extension of ISO_C_BINDING to be used if the
dnl#  gfortran extensions are available.
AC_DEFUN([WO_FC_ISO_C_BINDING_GFORTRAN_EMPTY],
 [module iso_c_binding_gfortran
  end module iso_c_binding_gfortran])

dnl#  Check for gfortran extensions in iso_c_binding
AC_DEFUN([WO_FC_CHECK_ISO_C_BINDING_GFORTRAN],
 [AC_CACHE_CHECK([whether $FC supports c_float128 (a gfortran extension)],
   [wo_cv_fc_iso_c_binding_gfortran],
   AC_REQUIRE([AC_PROG_FC])
   AC_LANG([Fortran])
   [AC_LINK_IFELSE(
      [dnl
       program conftest
         use iso_c_binding
         implicit none
         integer :: i
         i = c_float128
       end program conftest],
      [wo_cv_fc_iso_c_binding_gfortran=yes],
      [wo_cv_fc_iso_c_binding_gfortran=no])])])

dnl#  A, by autoconf standards, large program that performs
dnl#  RUNTIME checks of the available kinds.
dnl#  Note that it can NOT be used for cross compiling!
AC_DEFUN([WO_FC_MODULE_QUERY_KINDS],
[dnl
module query_kinds
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: iso_c_binding
  use :: iso_c_binding_gfortran
  implicit none
  private
  type, public :: real_kind
     integer :: kind = -1
     integer :: min_prec = -1
     integer :: max_prec = -1
     character(len=9) :: name  = "-1"
     character(len=7) :: iso_name = "-1"
     character(len=13) :: c_name = "-1"
     character(len=21) :: c_name_complex = "-1"
  end type real_kind
  type, public :: int_kind
     integer :: kind = -1
     integer :: min_range = -1
     integer :: max_range = -1
     character(len=9) :: name = "-1"
     character(len=5) :: iso_name = "-1"
  end type int_kind
  public :: query_real
  public :: query_int
  public :: real_by_name
  integer, parameter :: NUM_KINDS = 20
  integer, parameter :: MAX_PREC  = 99
  integer, parameter :: MAX_RANGE = 99
contains
  subroutine real_iso_name (kind)
    type(real_kind), intent(inout) :: kind
    !!! The following MUST NOT be a SELECT CASE statement,
    !!! because two or more of the real<n> kinds might be
    !!! negative and identical
    if (real32 >= 0 .and. kind%kind == real32) then
       kind%iso_name = "real32"
    else if (real64 >= 0 .and. kind%kind == real64) then
       kind%iso_name = "real64"
    else if (real128 >= 0 .and. kind%kind == real128) then
       kind%iso_name = "real128"
    end if
  end subroutine real_iso_name
  subroutine real_c_name (kind)
    type(real_kind), intent(inout) :: kind
    !!! The following MUST NOT be a SELECT CASE statement.
    if (c_float >= 0 .and. kind%kind == c_float) then
       kind%c_name = "c_float"
       kind%c_name_complex = "c_float_complex"
    else if (c_double >= 0 .and. kind%kind == c_double) then
       kind%c_name = "c_double"
       kind%c_name_complex = "c_double_complex"
    else if (c_long_double >= 0 .and. kind%kind == c_long_double) then
       kind%c_name = "c_long_double"
       kind%c_name_complex = "c_long_double_complex"
    !!! The following is gfortran specific:
    else if (c_float128 >= 0 .and. kind%kind == c_float128) then
       kind%c_name = "c_float128"
       kind%c_name_complex = "c_float128_complex"
    end if
  end subroutine real_c_name
  subroutine int_iso_name (kind)
    type(int_kind), intent(inout) :: kind
    !!! The following MUST NOT be a SELECT CASE statement.
    if (int8 >= 0 .and. kind%kind == int8) then
       kind%iso_name = "int8"
    else if (int16 >= 0 .and. kind%kind == int16) then
       kind%iso_name = "int16"
    else if (int32 >= 0 .and. kind%kind == int32) then
       kind%iso_name = "int32"
    else if (int64 >= 0 .and. kind%kind == int64) then
       kind%iso_name = "int64"
    end if
  end subroutine int_iso_name
  subroutine query_real (single, double, other)
    type(real_kind), intent(out) :: single
    type(real_kind), intent(out) :: double
    type(real_kind), dimension(:), allocatable, intent(out) :: other
    type(real_kind), dimension(NUM_KINDS) :: kinds
    integer :: kind_single = kind (1.0)
    integer :: kind_double = kind (1.0D0)
    integer :: precision_double = precision (1.0D0)
    integer :: p, k, last_k, offset, i, j
    single = real_kind ()
    double = real_kind ()
    last_k = -1
    offset = 0
    prec_loop: do p = 1, MAX_PREC
       k = selected_real_kind (p = p)
       if (k < 0) then
          exit prec_loop
       end if
       if (k == last_k) then
          kinds(offset)%max_prec = p
       else
          last_k = k
          offset = offset + 1
          kinds(offset)%kind = k
          kinds(offset)%min_prec = p
          kinds(offset)%max_prec = p
       end if
    end do prec_loop
    allocate (other(offset-2))
    i = 1
    do j = 1, offset
       if (kinds(j)%kind == kind_single) then
          single = kinds(j)
          single%name = "single"
          call real_iso_name (single)
          call real_c_name (single)
       else if (kinds(j)%kind == kind_double) then
          double = kinds(j)
          double%name = "double"
          call real_iso_name (double)
          call real_c_name (double)
       else
          if (i > offset - 2) then
             print *, "query_real: expected REAL and DOUBLE"
             return
          end if
          other(i) = kinds(j)
          if (other(i)%max_prec >= 4 * precision_double) then
             other(i)%name = "octuple"
          else if (other(i)%max_prec >= 2 * precision_double) then
             other(i)%name = "quadruple"
          else if (other(i)%max_prec > precision_double) then
             other(i)%name = "extended"
          else
             write (other(i)%name, "('prec',I2.2)") other(i)%max_prec
          endif
          call real_iso_name (other(i))
          call real_c_name (other(i))
          i = i + 1
       end if
    end do
  end subroutine query_real
  function name_matches (name, k) result (yorn)
    character(len=*), intent(in) :: name
    type(real_kind), intent(in) :: k
    logical :: yorn
    yorn = trim (name) == k%name &
         .or. trim (name) == k%iso_name &
         .or. trim (name) == k%c_name
  end function name_matches
  function real_by_name (name, single, double, other) result (match)
    character(len=*), intent(in) :: name
    type(real_kind), intent(in) :: single
    type(real_kind), intent(in) :: double
    type(real_kind), dimension(:), intent(in) :: other
    type(real_kind) :: match
    integer :: i
    match = real_kind ()
    if (name_matches (name, single)) then
       match = single
    else if (name_matches (name, double)) then
       match = double
    else
       do i = 1, size (other)
          if (name_matches (name, other(i))) then
             match = other(i)
             match%name = name
          end if
       end do
    end if
  end function real_by_name
  subroutine query_int (default, other)
    type(int_kind), intent(out) :: default
    type(int_kind), dimension(:), allocatable, intent(inout) :: other
    type(int_kind), dimension(NUM_KINDS) :: kinds
    integer :: kind_default = kind (1)
    integer :: r, k, last_k, offset, i, j
    default = int_kind ()
    last_k = -1
    offset = 0
    range_loop: do r = 1, MAX_PREC
       k = selected_int_kind (r)
       if (k < 0) then
          exit range_loop
       end if
       if (k == last_k) then
          kinds(offset)%max_range = r
       else
          last_k = k
          offset = offset + 1
          kinds(offset)%kind = k
          kinds(offset)%min_range = r
          kinds(offset)%max_range = r
       end if
       if (k == kind_default) then
          default = kinds(offset)
       end if
    end do range_loop
    allocate (other(offset-1))
    i = 1
    do j = 1, offset
       if (kinds(j)%kind == kind_default) then
          default = kinds(j)
          default%name = "dflt_int"
          call int_iso_name (default)
       else
          if (i > offset - 1) then
             print *, "query_int: expected DEFAULT"
             return
          end if
          other(i) = kinds(j)
          write (other(i)%name, "('range',I2.2)") other(i)%max_range
          call int_iso_name (other(i))
          i = i + 1
       end if
    end do
  end subroutine query_int
end module query_kinds])

AC_DEFUN([WO_FC_MODULE_REPORT_KINDS],
[dnl
module report_kinds
  use query_kinds
  implicit none
  private
  public print_kinds
contains
  subroutine print_kinds (request, single_real, double_real, other_real, &
                          default_int, other_int)
    character(len=*), intent(in)  :: request
    type(real_kind), intent(in) :: single_real, double_real
    type(real_kind), dimension(:), intent(in) :: other_real
    type(int_kind), intent(in) :: default_int
    type(int_kind), dimension(:), intent(in), allocatable :: other_int
    type(real_kind) :: default_real
    integer :: i
    print *, "module kinds"
    print *, " use, intrinsic :: iso_fortran_env"
    print *, " use, intrinsic :: iso_c_binding"
    print *, " implicit none"
    print *, " private"
    print *, ""
    print *, " !!! available REAL kinds               ! prec.  ! ISO     ! C"
    call print_real_kind (single_real)
    call print_real_kind (double_real)
    do i = 1, size (other_real)
       call print_real_kind (other_real(i))
    end do
    print *, ""
    print *, " !!! available INTEGER kinds            ! range  ! ISO     ! C"
    call print_int_kind (default_int)
    do i = 1, size (other_int)
       call print_int_kind (other_int(i))
    end do
    print *, ""
    print *, " !!! additional INTEGER kinds"
    print *, " public :: i8, i16, i32, i64"
    print *, " integer, parameter :: i8  = selected_int_kind (2)"
    print *, " integer, parameter :: i16 = selected_int_kind (4)"
    print *, " integer, parameter :: i32 = selected_int_kind (9)"
    print *, " integer, parameter :: i64 = selected_int_kind (18)"
    print *, " public :: TC"
    print *, " integer, parameter :: TC = i32"
    default_real = real_by_name (request, single_real, double_real, other_real)
    print *, ""
    print *, " !!! default REAL kinds"
    print *, " public :: single, double"
    print *, " public :: default, c_default_float, c_default_complex"
    print *, " integer, parameter :: default           = ", default_real%name
    print *, " integer, parameter :: c_default_float   = ", default_real%c_name
    print *, " integer, parameter :: c_default_complex = ", &
         default_real%c_name_complex
    print *, ""
    if (default_real%kind < 0) then
       print *, " !!! ERROR: the requested default real kind '" // &
               request // "' is not available in this environment"
       print *, ""
    else
       if (trim (default_real%c_name) == "-1" .or. &
         trim (default_real%c_name_complex) == "-1") then
       	   print *, " !!! ERROR: for the requested default real kind '" // &
           	     request // "' there is no supported C analogue"
           print *, ""
       else if (trim (request) == "real128") then
          if (default_real%max_prec < 2 * double_real%max_prec) then
             print *, " !!! WARNING: the requested default real kind 'real128'", &
               " does NOT provide quadruple precision in this environment"
             print *, ""
          end if
       end if
    end if
    print *, "end module kinds"
  end subroutine print_kinds
  subroutine print_real_kind (k)
    type(real_kind), intent(in) :: k
    write (*, "(A,A,A,I2,A,I2,A,I2,A,A,A,A)") &
         "  integer, parameter :: ", k%name, " = ", k%kind, &
         "   ! ", k%min_prec, "..", k%max_prec, &
         " ! ", k%iso_name, " ! ", k%c_name
  end subroutine print_real_kind
  subroutine print_int_kind (k)
    type(int_kind), intent(in) :: k
    write (*, "(A,A,A,I2,A,I2,A,I2,A,A,A,A)") &
         "  integer, parameter :: ", k%name, " = ", k%kind, &
         "   ! ", k%min_range, "..", k%max_range, &
         " ! ", k%iso_name, "   ! ", ""
  end subroutine print_int_kind
end module report_kinds])

AC_DEFUN([WO_FC_PROGRAM_CONFIGURE_KINDS],
[dnl
program configure_kinds
  use query_kinds
  use report_kinds
  implicit none
  type(real_kind), dimension(:), allocatable :: other_real
  type(int_kind), dimension(:), allocatable :: other_int
  type(real_kind) :: single_real, double_real
  type(int_kind) :: default_int
  character(len=120) :: request
  integer :: length, status
  call get_command_argument (1, request, length, status)
  if (status == 0 .and. length > 0) then
     call query_real (single_real, double_real, other_real)
     call query_int (default_int, other_int)
     call print_kinds (request(1:length), &
                       single_real, double_real, other_real, &
                       default_int, other_int)
  end if
end program configure_kinds])

dnl#  combining the source code
AC_DEFUN([WO_FC_CONFIGURE_KINDS_SOURCE],
[dnl
WO_FC_MODULE_QUERY_KINDS
WO_FC_MODULE_REPORT_KINDS
WO_FC_PROGRAM_CONFIGURE_KINDS])

AC_DEFUN([WO_FC_MODULE_KINDS_CROSS_COMPILING],
[dnl
! A minimalistic kinds.f90, since we're cross compiling
module kinds
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: single, double
  integer, parameter :: single = kind (1.0)
  integer, parameter :: double = kind (1.0D0)
  public :: i8, i16, i32, i64, TC
  integer, parameter :: i8  = selected_int_kind (2)
  integer, parameter :: i16 = selected_int_kind (4)
  integer, parameter :: i32 = selected_int_kind (9)
  integer, parameter :: i64 = selected_int_kind (18)
  integer, parameter :: TC = i32
  public :: default
  integer, parameter :: default = $1
  public :: c_default_float, c_default_complex
  integer, parameter :: c_default_float = $2
  integer, parameter :: c_default_complex = $3
end module kinds])

dnl#  fall back option for cross compilation
AC_DEFUN([WO_FC_CONFIGURE_KINDS_CROSS_COMPILING],
[dnl
  case "X$wo_cv_fc_requested_precision" in
    Xsingle)
      wo_cv_fc_requested_default_c_float_kind=c_float
      wo_cv_fc_requested_default_c_complex_kind=c_float_complex
      ;;
    Xdouble)
      wo_cv_fc_requested_default_c_float_kind=c_double
      wo_cv_fc_requested_default_c_complex_kind=c_double_complex
      ;;
    *)
      WO_FC_MSG_ERROR_BOX(dnl
        [cross compilation supported for single and double precision only!])
      ;;
  esac
  AC_COMPILE_IFELSE(dnl
    [WO_FC_MODULE_KINDS_CROSS_COMPILING(dnl
       [$wo_cv_fc_requested_precision],
       [$wo_cv_fc_requested_default_c_float_kind],
       [$wo_cv_fc_requested_default_c_complex_kind])],
    [cp conftest.$ac_ext kinds.f90],
    [WO_FC_MSG_ERROR_BOX(dnl
      [setting up kinds.f90 for cross compilation failed])])])

dnl#  report success
AC_DEFUN([WO_FC_CONFIGURE_KINDS_RUN_OK],
[dnl
  $INSTALL -d `AS_DIRNAME(["$1"])`
  ./conftest$EXEEXT "$wo_cv_fc_requested_precision" >$1
  if test X"$KEEP_KINDLY" != X; then
     cp -a ./conftest.f90 "$KEEP_KINDLY"
  fi
  if test ! -s $1; then
     WO_FC_MSG_ERROR_BOX([./conftest$EXEEXT produced no output])
  fi
  kinds_error="`$SED -n '/ERROR/s/^.*ERROR: *//p' $1`"
  if test "x$kinds_error" != x; then
     WO_FC_MSG_ERROR_BOX([$kinds_error])
  fi
  kinds_warning="`$SED -n '/WARNING/s/^.*WARNING: *//p' $1`"
  if test "x$kinds_warning" != x; then
     WO_FC_MSG_WARN_BOX([$kinds_warning])
  fi])

dnl#  report failure
AC_DEFUN([WO_FC_CONFIGURE_KINDS_RUN_FAIL],
[rm -f $1
 if test X"$KEEP_KINDLY" != X; then
    cp -a ./conftest.f90 "$KEEP_KINDLY"
 fi
 WO_FC_MSG_ERROR_BOX([could not compile kinds selection program])])

########################################################################
###
###  Prepare a kinds.f90 with a DEFAULT kind taken from the argument
###  of `--with-precision'.  If we are NOT cross compiling, 
###  add comments explaing the relationship to the kinds defined in
###  ISO_FORTRAN_ENV and ISO_C_BINDING.
###
###  The default of DEFAULT is `double'.
###
########################################################################
AC_DEFUN([WO_FC_CONFIGURE_KINDS],
[dnl
AC_REQUIRE([AC_PROG_SED])
AC_REQUIRE([AC_PROG_INSTALL])
AC_REQUIRE([WO_FC_FILENAME_CASE_CONVERSION])
AC_REQUIRE([AC_PROG_FC])
AC_LANG([Fortran])
AC_REQUIRE([WO_FC_CHECK_ISO_C_BINDING_GFORTRAN])
AC_MSG_CHECKING([the requested floating point precision])
wo_cv_fc_requested_precision=double
AC_ARG_WITH([precision],
  [  --with-precision=single|double|quadruple|extended|real32|real64|real128
                          request a floating point precision other than
                          double precision.  Note that only single and
                          double are guaranteed to be provided by all
                          Fortran compilers.],
  [case "x$withval" in
     x | xno | xyes ) wo_cv_fc_requested_precision=double ;;
     * )              wo_cv_fc_requested_precision="`echo $withval | $LOWERCASE`" ;;
   esac])
case "$wo_cv_fc_requested_precision" in
   single | double | extended | quadruple | real32 | real64 | real128 )
     AC_MSG_RESULT([$wo_cv_fc_requested_precision])
     ;;
   *)
     AC_MSG_RESULT()
     WO_FC_MSG_ERROR_BOX2([argument of --with-precision is $wo_cv_fc_requested_precision, but must be one of],
       [single, double, quadruple, extended, real32, real64, real128!])
     ;;
esac

dnl  save_cross_compiling=$cross_compiling
dnl  cross_compiling=yes

if test "x$wo_cv_fc_iso_c_binding_gfortran" = xyes; then
    AC_RUN_IFELSE(dnl
     [WO_FC_ISO_C_BINDING_GFORTRAN_EMPTY
      WO_FC_CONFIGURE_KINDS_SOURCE],
     [WO_FC_CONFIGURE_KINDS_RUN_OK([$1])],
     [WO_FC_CONFIGURE_KINDS_RUN_FAIL([$1])],
     [WO_FC_CONFIGURE_KINDS_CROSS_COMPILING])
else
    AC_RUN_IFELSE(dnl
     [WO_FC_ISO_C_BINDING_GFORTRAN_DUMMY
      WO_FC_CONFIGURE_KINDS_SOURCE],
     [WO_FC_CONFIGURE_KINDS_RUN_OK([$1])],
     [WO_FC_CONFIGURE_KINDS_RUN_FAIL([$1])],
     [WO_FC_CONFIGURE_KINDS_CROSS_COMPILING])
fi
rm -f iso_c_binding_gfortran.*
rm -f query_kinds.*
rm -f report_kinds.*

dnl  cross_compiling=$save_cross_compiling

if test "$wo_cv_fc_requested_precision" = "real128" -a "$FC_VENDOR" = "gfortran"; then
  FC_PRECISION="extended"
elif test "$wo_cv_fc_requested_precision" = "real128" -a "$FC_VENDOR" = "Intel"; then
  FC_PRECISION="quadruple"
elif test "$wo_cv_fc_requested_precision" = "real64"; then
  FC_PRECISION="double"
elif test "$wo_cv_fc_requested_precision" = "real32"; then
  FC_PRECISION="single"
else
  FC_PRECISION=$wo_cv_fc_requested_precision
fi
AC_SUBST([FC_PRECISION])
AM_CONDITIONAL([FC_PREC], [test "$FC_PRECISION" = "extended" || test "$FC_PRECISION" = "quadruple"])
AM_CONDITIONAL([FC_EXT], [test "$FC_PRECISION" = "extended"])
AM_CONDITIONAL([FC_QUAD], [test "$FC_PRECISION" = "quadruple"])
])

########################################################################
### end of configure kinds.f90
########################################################################


### filename_case_conversion, define two variables LOWERCASE and 
### UPPERCASE for /bin/sh filters that convert strings to lower 
### and upper case, respectively
AC_DEFUN([WO_FC_FILENAME_CASE_CONVERSION],
[dnl
AC_SUBST([LOWERCASE])
AC_SUBST([UPPERCASE])
AC_PATH_PROGS(TR,tr)
AC_MSG_CHECKING([for case conversion])
if test -n "$TR"; then
  LOWERCASE="$TR A-Z a-z"
  UPPERCASE="$TR a-z A-Z"
  WO_FC_FILENAME_CASE_CONVERSION_TEST
fi
if test -n "$UPPERCASE" && test -n "$LOWERCASE"; then
  AC_MSG_RESULT([$TR works])
else
  LOWERCASE="$SED y/ABCDEFGHIJKLMNOPQRSTUVWXYZ/abcdefghijklmnopqrstuvwxyz/"
  UPPERCASE="$SED y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/"
  WO_FC_FILENAME_CASE_CONVERSION_TEST
  if test -n "$UPPERCASE" && test -n "$LOWERCASE"; then
    AC_MSG_RESULT([$SED works])
  fi
fi])
### end WO_FC_FILE_CASE_CONVERSION
dnl
AC_DEFUN([WO_FC_FILENAME_CASE_CONVERSION_TEST],
[dnl
if test "`echo fOo | $LOWERCASE`" != "foo"; then
  LOWERCASE=""
fi
if test "`echo fOo | $UPPERCASE`" != "FOO"; then
  UPPERCASE=""
fi])

dnl
dnl --------------------------------------------------------------------
dnl
dnl COMPILE_FC(VARIABLE, COMPILER, EXTENSION, MODULE,
dnl                       VALUE_SUCCESS, VALUE_FAILURE, KEEP)
dnl
AC_DEFUN([COMPILE_FC],
[cat >conftest.$3 <<__END__
$4
program conftest
  print *, 42
end program conftest
__END__
$2 -o conftest conftest.$3 >/dev/null 2>&1
./conftest >conftest.out 2>/dev/null
if test 42 = "`$SED 's/ //g' conftest.out`"; then
  $1="$5"
else
  $1="$6"
fi
if test -z "$7"; then
  rm -rf conftest* CONFTEST*
fi])

dnl
dnl --------------------------------------------------------------------
dnl
dnl WO_FC_MODULE_FILE(NAME, EXTENSION, FORTRAN_COMPILER, SOURCE_EXTENSION)
dnl
AC_DEFUN([WO_FC_MODULE_FILE],
[AC_SUBST([$1])
AC_SUBST([$2])
AC_MSG_CHECKING([for Fortran90 module file naming convention])
COMPILE_FC([wo_result], [$3], [$4],
  [module module_NAME
     implicit none
     integer, parameter, public :: forty_two = 42
   end module module_NAME], [ok], [], [KEEP])
if test -n "[$]wo_result"; then
  $1=unknown
  $2=unknown
  for name in module_name module_NAME MODULE_NAME conftest; do
    for ext in m mod M MOD d D; do
      if test -f "[$]name.[$]ext"; then
        $1="$name"
        $2="$ext"
        break 2
      fi
    done
  done
  if test X"[$]$1" = X"module_name"; then
     AC_MSG_RESULT([name: [$]$1, extension: .[$]$2 ])
  else
     AC_MSG_ERROR([unusual unsupported module file name convention: [$]$1.[$]$2])
  fi
else
  $1=""
  $2=""
  AC_MSG_RESULT([compiler failed])
  AC_MSG_NOTICE([error: **************************************************************])
  AC_MSG_NOTICE([error: Fortran compiler cannot create proper module files. This])
  AC_MSG_NOTICE([error: might be caused by linking against a wrong gcc library.])
  AC_MSG_ERROR([**************************************************************])
fi
rm -rf conftest* CONFTEST* module_name* module_NAME* MODULE_NAME*])
